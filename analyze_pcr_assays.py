#!/usr/bin/env python3
"""
Script to analyze PCR assay designs and generate statistics to assess their validity.
"""

import re
import csv
import os
import primer3  # For accurate Tm calculations matching original design
from Bio.SeqUtils import gc_fraction  # Using Bio.SeqUtils.gc_fraction as per memory

def parse_block(block, block_index=0):
    """Parse a single PCR assay block into a structured dictionary."""
    # Skip empty blocks and the file header
    if not block.strip() or block.strip().startswith('# PCR Primer Design Results'):
        return None
    
    # Check if this block is a valid entry (must have forward & reverse primers)
    if not ('Forward Primer:' in block and 'Reverse Primer:' in block):
        return None
    
    # Create a new assay entry
    assay = {}
    
    # Parse forward primer
    fw_match = re.search(r'Forward Primer: (\w+) \(Tm: ([\d.]+), GC: ([\d.]+)%\)', block)
    if fw_match:
        assay['forward_primer'] = fw_match.group(1)
        assay['forward_tm'] = float(fw_match.group(2))
        assay['forward_gc'] = float(fw_match.group(3))
    else:
        return None  # Skip if we can't parse forward primer
    
    # Parse reverse primer
    rev_match = re.search(r'Reverse Primer: (\w+) \(Tm: ([\d.]+), GC: ([\d.]+)%\)', block)
    if rev_match:
        assay['reverse_primer'] = rev_match.group(1)
        assay['reverse_tm'] = float(rev_match.group(2))
        assay['reverse_gc'] = float(rev_match.group(3))
    else:
        return None  # Skip if we can't parse reverse primer
    
    # Parse amplicon sequence and length
    amp_match = re.search(r'Amplicon Sequence: (\w+) \(Length: (\d+)\)', block)
    if amp_match:
        assay['amplicon'] = amp_match.group(1)
        assay['amplicon_length'] = int(amp_match.group(2))
    
    # Parse amplicon reverse complement
    amp_rc_match = re.search(r'Amplicon Reverse Complement: (\w+)', block)
    if amp_rc_match:
        assay['amplicon_rc'] = amp_rc_match.group(1)
    
    # Parse original target
    target_match = re.search(r'Original Target: (\w+)', block)
    if target_match:
        assay['original_target'] = target_match.group(1)
    
    # Parse metadata and extract virus information
    meta_match = re.search(r'Metadata: ([^\n]+)', block)
    if meta_match:
        metadata = meta_match.group(1)
        assay['metadata'] = metadata
        
        # Extract virus name from metadata
        virus_match = re.search(r'\|(.*?)\|', metadata)
        if virus_match:
            assay['virus'] = virus_match.group(1).strip()
    
    return assay

def process_chunk(chunk_blocks, chunk_id):
    """Process a chunk of blocks in parallel."""
    chunk_assays = []
    skipped = 0
    
    for i, block in enumerate(chunk_blocks):
        assay = parse_block(block, i)
        if assay:
            chunk_assays.append(assay)
        else:
            skipped += 1
    
    return {'assays': chunk_assays, 'skipped': skipped}

def parse_pcr_assays(file_path):
    """Parse PCR assay designs from the input file using multiprocessing."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    
    # Check file size first to avoid memory issues
    file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB
    print(f"File size: {file_size:.2f} MB")
    
    # Read the entire file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Split on double newlines to separate entries
    entry_blocks = re.split(r'\n\s*\n', content.strip())
    print(f"Found {len(entry_blocks)} blocks in the file")
    
    # Determine optimal chunk size and number of processes
    num_cpus = mp.cpu_count()
    num_processes = max(1, num_cpus - 1)  # Leave one CPU for system tasks
    
    total_blocks = len(entry_blocks)
    chunk_size = max(1, total_blocks // (num_processes * 10))  # Create ~10 chunks per process
    
    print(f"Using {num_processes} processes with chunk size of {chunk_size}")
    
    # Show sample entries
    for i in range(min(3, len(entry_blocks))):
        sample_assay = parse_block(entry_blocks[i], i)
        if sample_assay:
            print(f"Sample entry {i + 1}:")
            print(f"  Forward primer: {sample_assay.get('forward_primer', 'N/A')}")
            print(f"  Reverse primer: {sample_assay.get('reverse_primer', 'N/A')}")
            print(f"  Amplicon length: {sample_assay.get('amplicon_length', 'N/A')}")
    
    # Create chunks of blocks to process
    chunks = [entry_blocks[i:i+chunk_size] for i in range(0, len(entry_blocks), chunk_size)]
    print(f"Created {len(chunks)} chunks for processing")
    
    # Process chunks in parallel
    assays = []
    skipped_blocks = 0
    
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = [executor.submit(process_chunk, chunk, i) for i, chunk in enumerate(chunks)]
        
        for future in futures:
            result = future.result()
            assays.extend(result['assays'])
            skipped_blocks += result['skipped']
    
    print(f"Skipped {skipped_blocks} blocks out of {len(entry_blocks)}")
    print(f"Successfully parsed {len(assays)} PCR assay designs")
    
    return assays

def calculate_statistics(assays):
    """Calculate statistics for PCR assay designs."""
    import numpy as np
    from collections import defaultdict
    
    # Process each assay to add derived metrics
    for assay in assays:
        # Calculate derived metrics
        if 'forward_primer' in assay and 'reverse_primer' in assay:
            assay['forward_length'] = len(assay['forward_primer'])
            assay['reverse_length'] = len(assay['reverse_primer'])
            assay['tm_difference'] = abs(assay['forward_tm'] - assay['reverse_tm'])
            
            # Check for primer dimers (simplified check for complementary sequences at 3' end)
            primer_3prime_fw = assay['forward_primer'][-5:]
            primer_3prime_rev = assay['reverse_primer'][-5:]
            rev_comp_fw = ''.join([{'A':'T', 'T':'A', 'G':'C', 'C':'G'}.get(base, base) for base in primer_3prime_fw[::-1]])
            rev_comp_rev = ''.join([{'A':'T', 'T':'A', 'G':'C', 'C':'G'}.get(base, base) for base in primer_3prime_rev[::-1]])
            
            # Check if at least 3 bases at 3' end are complementary
            assay['potential_self_dimer_fw'] = primer_3prime_fw[-3:] in rev_comp_fw
            assay['potential_self_dimer_rev'] = primer_3prime_rev[-3:] in rev_comp_rev
            assay['potential_cross_dimer'] = primer_3prime_fw[-3:] in rev_comp_rev or primer_3prime_rev[-3:] in rev_comp_fw
            
            # Verify GC content using Biopython's gc_fraction (multiplied by 100 for percentage)
            assay['forward_gc_calculated'] = gc_fraction(assay['forward_primer']) * 100
            assay['reverse_gc_calculated'] = gc_fraction(assay['reverse_primer']) * 100
            
            # Calculate melting temperatures using primer3's calc_tm function with default parameters
            try:
                # Using primer3.calc_tm instead of deprecated calcTm function
                # Using default parameters without specifying any arguments
                assay['forward_tm_calculated'] = primer3.calc_tm(assay['forward_primer'])
                assay['reverse_tm_calculated'] = primer3.calc_tm(assay['reverse_primer'])
                
                # Check if the provided Tm values are reasonably close to calculated values
                # Allow for some variance (±2°C) due to potentially different calculation methods
                assay['forward_tm_accurate'] = abs(assay['forward_tm'] - assay['forward_tm_calculated']) < 2.0
                assay['reverse_tm_accurate'] = abs(assay['reverse_tm'] - assay['reverse_tm_calculated']) < 2.0
                assay['tm_values_accurate'] = assay['forward_tm_accurate'] and assay['reverse_tm_accurate']
            except Exception as e:
                # Graceful handling of any calculation issues
                assay['forward_tm_calculated'] = 0.0
                assay['reverse_tm_calculated'] = 0.0
                assay['tm_values_accurate'] = False
            
        if 'amplicon' in assay and 'original_target' in assay:
            # Check if amplicon is contained within the original target
            assay['amplicon_in_target'] = assay['amplicon'] in assay['original_target']
            
            # Calculate amplicon percentage of the original target
            if 'original_target' in assay:
                assay['amplicon_percentage'] = (assay['amplicon_length'] / len(assay['original_target'])) * 100
        
        # Evaluate design quality
        assay['valid_tm_diff'] = assay.get('tm_difference', float('inf')) <= 5.0
        assay['valid_gc_content'] = (40 <= assay.get('forward_gc', 0) <= 60 and 
                                    40 <= assay.get('reverse_gc', 0) <= 60)
        assay['valid_primer_length'] = (18 <= assay.get('forward_length', 0) <= 25 and 
                                       18 <= assay.get('reverse_length', 0) <= 25)
        assay['valid_amplicon_size'] = 80 <= assay.get('amplicon_length', 0) <= 150
        assay['valid_tm_range'] = (55 <= assay.get('forward_tm', 0) <= 65 and
                                  55 <= assay.get('reverse_tm', 0) <= 65)
        assay['no_primer_dimers'] = not (assay.get('potential_self_dimer_fw', True) or 
                                         assay.get('potential_self_dimer_rev', True) or 
                                         assay.get('potential_cross_dimer', True))
        
        # Add Tm accuracy to validity assessment
        assay['valid_tm_accuracy'] = assay.get('tm_values_accurate', False)
        
        # Overall validity score
        validity_factors = [
            assay.get('valid_tm_diff', False),
            assay.get('valid_gc_content', False),
            assay.get('valid_primer_length', False),
            assay.get('valid_amplicon_size', False),
            assay.get('amplicon_in_target', False),
            assay.get('valid_tm_range', False),
            assay.get('no_primer_dimers', False),
            assay.get('valid_tm_accuracy', False)
        ]
        assay['validity_score'] = sum(validity_factors) / len(validity_factors) * 100
    
    # Calculate aggregate statistics after processing all assays
    stats = {}
    
    # Valid designs count
    valid_assays = [a for a in assays if a.get('validity_score', 0) >= 80]
    stats['total_valid_designs'] = len(valid_assays)
    
    # Extract arrays for numerical properties
    amplicon_lengths = np.array([a.get('amplicon_length', 0) for a in assays if 'amplicon_length' in a])
    fw_lengths = np.array([a.get('forward_length', 0) for a in assays if 'forward_length' in a])
    fw_gc = np.array([a.get('forward_gc', 0) for a in assays if 'forward_gc' in a])
    fw_tm = np.array([a.get('forward_tm', 0) for a in assays if 'forward_tm' in a])
    rev_lengths = np.array([a.get('reverse_length', 0) for a in assays if 'reverse_length' in a])
    rev_gc = np.array([a.get('reverse_gc', 0) for a in assays if 'reverse_gc' in a])
    rev_tm = np.array([a.get('reverse_tm', 0) for a in assays if 'reverse_tm' in a])
    
    # Calculate additional statistics for Tm calculations
    fw_tm_calculated = np.array([a.get('forward_tm_calculated', 0) for a in assays if 'forward_tm_calculated' in a and a['forward_tm_calculated'] > 0])
    rev_tm_calculated = np.array([a.get('reverse_tm_calculated', 0) for a in assays if 'reverse_tm_calculated' in a and a['reverse_tm_calculated'] > 0])
    tm_accurate_count = sum(1 for a in assays if a.get('tm_values_accurate', False))
    
    # Calculate statistics if we have data
    if len(amplicon_lengths) > 0:
        stats['amplicon_length_avg'] = np.mean(amplicon_lengths)
        stats['amplicon_length_std'] = np.std(amplicon_lengths)
        stats['amplicon_length_min'] = np.min(amplicon_lengths)
        stats['amplicon_length_max'] = np.max(amplicon_lengths)
        
    if len(fw_tm_calculated) > 0:
        stats['forward_tm_calculated_avg'] = np.mean(fw_tm_calculated)
        stats['forward_tm_calculated_std'] = np.std(fw_tm_calculated)
        stats['forward_tm_calculated_min'] = np.min(fw_tm_calculated)
        stats['forward_tm_calculated_max'] = np.max(fw_tm_calculated)
        stats['forward_tm_difference_avg'] = np.mean(np.abs(fw_tm - fw_tm_calculated))
        
    if len(rev_tm_calculated) > 0:
        stats['reverse_tm_calculated_avg'] = np.mean(rev_tm_calculated)
        stats['reverse_tm_calculated_std'] = np.std(rev_tm_calculated)
        stats['reverse_tm_calculated_min'] = np.min(rev_tm_calculated)
        stats['reverse_tm_calculated_max'] = np.max(rev_tm_calculated)
        stats['reverse_tm_difference_avg'] = np.mean(np.abs(rev_tm - rev_tm_calculated))
        
    stats['tm_accurate_count'] = tm_accurate_count
    stats['tm_accurate_percent'] = (tm_accurate_count / len(assays)) * 100 if len(assays) > 0 else 0
    
    if len(fw_lengths) > 0:
        stats['forward_length_avg'] = np.mean(fw_lengths)
        stats['forward_length_std'] = np.std(fw_lengths)
        stats['forward_length_min'] = np.min(fw_lengths)
        stats['forward_length_max'] = np.max(fw_lengths)
    
    if len(fw_gc) > 0:
        stats['forward_gc_avg'] = np.mean(fw_gc)
        stats['forward_gc_std'] = np.std(fw_gc)
        stats['forward_gc_min'] = np.min(fw_gc)
        stats['forward_gc_max'] = np.max(fw_gc)
    
    if len(fw_tm) > 0:
        stats['forward_tm_avg'] = np.mean(fw_tm)
        stats['forward_tm_std'] = np.std(fw_tm)
        stats['forward_tm_min'] = np.min(fw_tm)
        stats['forward_tm_max'] = np.max(fw_tm)
    
    if len(rev_lengths) > 0:
        stats['reverse_length_avg'] = np.mean(rev_lengths)
        stats['reverse_length_std'] = np.std(rev_lengths)
        stats['reverse_length_min'] = np.min(rev_lengths)
        stats['reverse_length_max'] = np.max(rev_lengths)
    
    if len(rev_gc) > 0:
        stats['reverse_gc_avg'] = np.mean(rev_gc)
        stats['reverse_gc_std'] = np.std(rev_gc)
        stats['reverse_gc_min'] = np.min(rev_gc)
        stats['reverse_gc_max'] = np.max(rev_gc)
    
    if len(rev_tm) > 0:
        stats['reverse_tm_avg'] = np.mean(rev_tm)
        stats['reverse_tm_std'] = np.std(rev_tm)
        stats['reverse_tm_min'] = np.min(rev_tm)
        stats['reverse_tm_max'] = np.max(rev_tm)
    
    return assays, stats

def save_to_csv(assays, output_file):
    """Save PCR assay statistics to a CSV file."""
    # Define fields for the CSV
    fieldnames = [
        'forward_primer', 'forward_length', 'forward_tm', 'forward_tm_calculated', 'forward_tm_accurate', 'forward_gc', 'forward_gc_calculated',
        'reverse_primer', 'reverse_length', 'reverse_tm', 'reverse_tm_calculated', 'reverse_tm_accurate', 'reverse_gc', 'reverse_gc_calculated',
        'tm_difference', 'valid_tm_diff', 'valid_tm_range', 'tm_values_accurate', 'valid_tm_accuracy',
        'amplicon_length', 'amplicon_in_target', 'amplicon_percentage',
        'valid_gc_content', 'valid_primer_length', 'valid_amplicon_size',
        'potential_self_dimer_fw', 'potential_self_dimer_rev', 'potential_cross_dimer', 'no_primer_dimers',
        'validity_score', 'virus'
    ]
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for assay in assays:
            # Create a row with only the fields we want in the CSV
            row = {field: assay.get(field, '') for field in fieldnames}
            writer.writerow(row)

def save_detailed_stats(stats, output_file):
    """Save detailed statistics to a CSV file."""
    import csv
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Statistic', 'Value'])
        
        # Write each statistic to the CSV
        for key, value in sorted(stats.items()):
            writer.writerow([key, f"{value:.4f}" if isinstance(value, float) else value])

def main():
    import time
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    
    start_time = time.time()
    
    input_file = "pcr_assays_output.txt"
    output_file = "pcr_assays_stats.csv"
    summary_file = "pcr_assays_summary.txt"
    detailed_stats_file = "pcr_assays_detailed_stats.csv"
    
    # Default to processing all entries, but can be limited for testing
    # None means process all entries
    max_entries = None
    
    # Performance optimization: determine number of processes based on CPU count
    num_processes = max(1, mp.cpu_count() - 1)  # Leave one CPU for system
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return
    
    # Parse PCR assays from input file
    print(f"Parsing PCR assays from '{input_file}'...")
    assays = parse_pcr_assays(input_file)
    parse_time = time.time()
    print(f"Parsing completed in {parse_time - start_time:.2f} seconds")
    
    # Limit the number of entries if specified
    if max_entries and len(assays) > max_entries:
        print(f"Limiting analysis to {max_entries} entries for testing purposes")
        assays = assays[:max_entries]
    
    if not assays:
        print("No PCR assays found in the input file.")
        return
    
    # Calculate statistics
    print(f"Calculating statistics for {len(assays)} entries...")
    assays_with_stats, stats = calculate_statistics(assays)
    calc_time = time.time()
    print(f"Statistics calculation completed in {calc_time - parse_time:.2f} seconds")
    
    # Save individual assay statistics to CSV
    print(f"Saving assay statistics to '{output_file}'...")
    save_to_csv(assays_with_stats, output_file)
    
    # Save detailed aggregate statistics
    print(f"Saving detailed statistics to '{detailed_stats_file}'...")
    save_detailed_stats(stats, detailed_stats_file)
    
    print(f"Analysis complete. Found {len(assays)} PCR assay designs.")
    print(f"Total valid designs: {stats['total_valid_designs']} ({stats['total_valid_designs']/len(assays)*100:.2f}%)")
    
    # Write summary to file
    with open(summary_file, 'w') as f:
        f.write(f"PCR Assay Design Statistics Summary\n")
        f.write(f"=================================\n\n")
        f.write(f"Total number of assay designs analyzed: {len(assays)}\n")
        f.write(f"Total valid designs: {stats['total_valid_designs']} ({stats['total_valid_designs']/len(assays)*100:.2f}%)\n\n")
        
        f.write(f"Amplicon Statistics:\n")
        f.write(f"  Average length: {stats['amplicon_length_avg']:.2f} bp\n")
        f.write(f"  Standard deviation: {stats['amplicon_length_std']:.2f} bp\n")
        f.write(f"  Minimum length: {stats['amplicon_length_min']} bp\n")
        f.write(f"  Maximum length: {stats['amplicon_length_max']} bp\n\n")
        
        f.write(f"Forward Primer Statistics:\n")
        f.write(f"  Average length: {stats['forward_length_avg']:.2f} bp\n")
        f.write(f"  Length standard deviation: {stats['forward_length_std']:.2f} bp\n")
        f.write(f"  Average GC content: {stats['forward_gc_avg']:.2f}%\n")
        f.write(f"  GC content standard deviation: {stats['forward_gc_std']:.2f}%\n")
        f.write(f"  Maximum GC content: {stats['forward_gc_max']:.2f}%\n")
        f.write(f"  Minimum GC content: {stats['forward_gc_min']:.2f}%\n")
        f.write(f"  Average melting temperature: {stats['forward_tm_avg']:.2f}°C\n")
        f.write(f"  Melting temperature standard deviation: {stats['forward_tm_std']:.2f}°C\n")
        f.write(f"  Maximum melting temperature: {stats['forward_tm_max']:.2f}°C\n")
        f.write(f"  Minimum melting temperature: {stats['forward_tm_min']:.2f}°C\n")
        if 'forward_tm_calculated_avg' in stats:
            f.write(f"  Average calculated Tm (Nearest-Neighbor): {stats['forward_tm_calculated_avg']:.2f}°C\n")
            f.write(f"  Average difference between provided and calculated Tm: {stats['forward_tm_difference_avg']:.2f}°C\n\n")
        else:
            f.write("\n")
        
        f.write(f"Reverse Primer Statistics:\n")
        f.write(f"  Average length: {stats['reverse_length_avg']:.2f} bp\n")
        f.write(f"  Length standard deviation: {stats['reverse_length_std']:.2f} bp\n")
        f.write(f"  Average GC content: {stats['reverse_gc_avg']:.2f}%\n")
        f.write(f"  GC content standard deviation: {stats['reverse_gc_std']:.2f}%\n")
        f.write(f"  Maximum GC content: {stats['reverse_gc_max']:.2f}%\n")
        f.write(f"  Minimum GC content: {stats['reverse_gc_min']:.2f}%\n")
        f.write(f"  Average melting temperature: {stats['reverse_tm_avg']:.2f}°C\n")
        f.write(f"  Melting temperature standard deviation: {stats['reverse_tm_std']:.2f}°C\n")
        f.write(f"  Maximum melting temperature: {stats['reverse_tm_max']:.2f}°C\n")
        f.write(f"  Minimum melting temperature: {stats['reverse_tm_min']:.2f}°C\n")
        if 'reverse_tm_calculated_avg' in stats:
            f.write(f"  Average calculated Tm (Nearest-Neighbor): {stats['reverse_tm_calculated_avg']:.2f}°C\n")
            f.write(f"  Average difference between provided and calculated Tm: {stats['reverse_tm_difference_avg']:.2f}°C\n")
        
        f.write(f"\nTm Accuracy:\n")
        f.write(f"  Primer pairs with accurate Tm values: {stats.get('tm_accurate_count', 0)} ({stats.get('tm_accurate_percent', 0):.2f}%)\n")
    
    end_time = time.time()
    print(f"Detailed summary saved to '{summary_file}'")
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
    print(f"Detailed statistics saved to '{detailed_stats_file}'")
    
    # Print the key statistics to console
    print("\nKey Statistics:")
    print(f"  Amplicon length:        {stats['amplicon_length_avg']:.2f} ± {stats['amplicon_length_std']:.2f} bp")
    print(f"  Forward primer length:  {stats['forward_length_avg']:.2f} ± {stats['forward_length_std']:.2f} bp")
    print(f"  Forward primer GC:      {stats['forward_gc_avg']:.2f} ± {stats['forward_gc_std']:.2f}%  (Min: {stats['forward_gc_min']:.2f}%, Max: {stats['forward_gc_max']:.2f}%)")
    print(f"  Forward primer Tm:      {stats['forward_tm_avg']:.2f} ± {stats['forward_tm_std']:.2f}°C  (Min: {stats['forward_tm_min']:.2f}°C, Max: {stats['forward_tm_max']:.2f}°C)")
    if 'forward_tm_calculated_avg' in stats:
        print(f"  Forward primer Tm (NN): {stats['forward_tm_calculated_avg']:.2f} ± {stats['forward_tm_calculated_std']:.2f}°C")
    print(f"  Reverse primer length:  {stats['reverse_length_avg']:.2f} ± {stats['reverse_length_std']:.2f} bp")
    print(f"  Reverse primer GC:      {stats['reverse_gc_avg']:.2f} ± {stats['reverse_gc_std']:.2f}%  (Min: {stats['reverse_gc_min']:.2f}%, Max: {stats['reverse_gc_max']:.2f}%)")
    print(f"  Reverse primer Tm:      {stats['reverse_tm_avg']:.2f} ± {stats['reverse_tm_std']:.2f}°C  (Min: {stats['reverse_tm_min']:.2f}°C, Max: {stats['reverse_tm_max']:.2f}°C)")
    if 'reverse_tm_calculated_avg' in stats:
        print(f"  Reverse primer Tm (NN): {stats['reverse_tm_calculated_avg']:.2f} ± {stats['reverse_tm_calculated_std']:.2f}°C")
    if 'tm_accurate_percent' in stats:
        print(f"  Tm accuracy rate:       {stats['tm_accurate_percent']:.2f}%")

if __name__ == "__main__":
    main()
