#!/usr/bin/env python3
"""
Analyze QDB assay designs from output file and generate statistical tables.
"""

import pandas as pd
import numpy as np
import re
import os
import csv
import glob
from collections import Counter
import matplotlib.pyplot as plt
from tabulate import tabulate

def extract_virus_genomes_from_fasta_directory(directory_path):
    """
    Extract virus genome information from FASTA files in the specified directory.
    
    Args:
        directory_path (str): Path to the directory containing aligned genome FASTA files
        
    Returns:
        dict: Dictionary mapping GenBank IDs and virus names to their corresponding file paths
    """
    if not os.path.exists(directory_path):
        print(f"WARNING: Directory not found: {directory_path}")
        return {}, {}
    
    genbank_id_map = {}  # Maps GenBank IDs to file paths
    virus_name_map = {}  # Maps virus names to file paths
    
    # Find all FASTA files in the directory
    fasta_files = glob.glob(os.path.join(directory_path, "*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files in {directory_path}")
    
    for fasta_file in fasta_files:
        try:
            with open(fasta_file, 'r') as f:
                header_line = f.readline().strip()
                
                # Header format: >GenBank_ID |Description|Virus_Name|...
                if header_line.startswith('>'):
                    header_parts = header_line[1:].split('|')
                    if len(header_parts) >= 3:
                        genbank_id = header_parts[0].strip()
                        virus_name = header_parts[2].strip()
                        
                        # Store mappings
                        genbank_id_map[genbank_id] = fasta_file
                        virus_name_map[virus_name] = fasta_file
        except Exception as e:
            print(f"Error processing file {fasta_file}: {e}")
    
    print(f"Extracted {len(genbank_id_map)} GenBank IDs and {len(virus_name_map)} virus names from FASTA files")
    return genbank_id_map, virus_name_map

def parse_qdb_design_file(file_path):
    """
    Parse the QDB assay output file and extract relevant information.
    
    Args:
        file_path (str): Path to the QDB assay output file
        
    Returns:
        list: List of dictionaries containing design information
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    designs = []
    current_design = {}
    valid_designs_count = 0
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if not line and current_design:  # Empty line and we have a design to add
            # Only add designs with complete info
            if ('target_sequence' in current_design and 
                'capture_probe' in current_design and 
                'reporter_probe' in current_design):
                valid_designs_count += 1
            designs.append(current_design)
            current_design = {}
        elif line.startswith("Target:"):
            current_design['target_sequence'] = line.replace("Target:", "").strip()
        elif line.startswith("Capture Probe:"):
            # Try multiple regex patterns to handle different formats
            match = re.search(r'/5AmMC6/(.+) \(GC: (.+)%, Tm: (.+)°C\)', line)
            if not match:
                # Try alternative pattern
                match = re.search(r'Capture Probe: (.+) \(GC: (.+)%, Tm: (.+)°C\)', line)
            
            if match:
                current_design['capture_probe'] = match.group(1)
                current_design['capture_gc'] = float(match.group(2))
                current_design['capture_tm'] = float(match.group(3))
        elif line.startswith("Reporter Probe:"):
            # Extract only the nucleotide sequence, excluding markers
            # First try to match the format with /3Cy5Sp/ marker
            match = re.search(r'Reporter Probe: ([ATGC]+)/3Cy5Sp/ \(GC: (.+)%, Tm: (.+)°C\)', line)
            if not match:
                # Try alternative format without marker
                match = re.search(r'Reporter Probe: ([ATGC]+) \(GC: (.+)%, Tm: (.+)°C\)', line)
                
            if match:
                current_design['reporter_probe'] = match.group(1)  # Only nucleotides
                current_design['reporter_gc'] = float(match.group(2))
                current_design['reporter_tm'] = float(match.group(3))
        elif line.startswith("GenBank ID:"):
            current_design['genbank_id'] = line.replace("GenBank ID:", "").strip()
        elif line.startswith("Metadata:"):
            metadata = line.replace("Metadata:", "").strip()
            current_design['metadata'] = metadata
            
            # Parse metadata for additional information
            metadata_parts = metadata.split("|")
            
            # Extract organism info - this typically contains strain information
            if len(metadata_parts) > 1:
                full_organism = metadata_parts[1].strip()
                current_design['full_organism_info'] = full_organism
                
                # Extract base organism name without strain
                organism_match = re.match(r'^([^,]+)(?:\s+strain\s+|\s+isolate\s+|\s*,)', full_organism)
                if organism_match:
                    current_design['organism'] = organism_match.group(1).strip()
                else:
                    current_design['organism'] = full_organism
            else:
                current_design['organism'] = "Unknown"
                current_design['full_organism_info'] = "Unknown"
            
            # Extract virus name (standardized name)
            if len(metadata_parts) > 2:
                current_design['virus_name'] = metadata_parts[2].strip()
            else:
                current_design['virus_name'] = "Unknown"
            
            # Extract proper taxonomy info
            if len(metadata_parts) > 8 and metadata_parts[8].strip():
                current_design['virus_genus'] = metadata_parts[8].strip()
            else:
                current_design['virus_genus'] = "Unknown"
                
            if len(metadata_parts) > 9 and metadata_parts[9].strip():
                current_design['virus_family'] = metadata_parts[9].strip()
            else:
                current_design['virus_family'] = "Unknown"
                
            if len(metadata_parts) > 10 and metadata_parts[10].strip():
                current_design['genome_type'] = metadata_parts[10].strip()
            else:
                current_design['genome_type'] = "Unknown"
        
        i += 1
        
    # Add the last design if it exists
    if current_design:
        # Only count designs with complete info
        if ('target_sequence' in current_design and 
            'capture_probe' in current_design and 
            'reporter_probe' in current_design):
            valid_designs_count += 1
        designs.append(current_design)
    
    print(f"Total designs parsed: {len(designs)}")
    print(f"Valid designs with complete probe information: {valid_designs_count}")
        
    return designs

def calculate_statistics(designs, genbank_id_map=None, virus_name_map=None):
    """
    Calculate statistics for the designs.
    
    Args:
        designs (list): List of dictionaries containing design information
        genbank_id_map (dict, optional): Dictionary mapping GenBank IDs to FASTA file paths
        virus_name_map (dict, optional): Dictionary mapping virus names to FASTA file paths
        
    Returns:
        dict: Dictionary containing calculated statistics
    """
    stats = {}
    
    # Filter only valid designs with complete probe information
    valid_designs = [d for d in designs if (
        'target_sequence' in d and 
        'capture_probe' in d and 
        'reporter_probe' in d and
        'capture_gc' in d and
        'reporter_gc' in d and
        'capture_tm' in d and
        'reporter_tm' in d
    )]
    
    if not valid_designs:
        print("ERROR: No valid designs found with complete probe information!")
        # Return default stats to avoid errors
        return {
            'total_designs': 0,
            'avg_capture_gc': 0, 'std_capture_gc': 0, 'min_capture_gc': 0, 'max_capture_gc': 0,
            'avg_reporter_gc': 0, 'std_reporter_gc': 0, 'min_reporter_gc': 0, 'max_reporter_gc': 0,
            'avg_capture_tm': 0, 'std_capture_tm': 0, 'min_capture_tm': 0, 'max_capture_tm': 0,
            'avg_reporter_tm': 0, 'std_reporter_tm': 0, 'min_reporter_tm': 0, 'max_reporter_tm': 0,
            'avg_target_length': 0, 'std_target_length': 0, 'min_target_length': 0, 'max_target_length': 0,
            'avg_capture_length': 0, 'std_capture_length': 0, 'min_capture_length': 0, 'max_capture_length': 0,
            'avg_reporter_length': 0, 'std_reporter_length': 0, 'min_reporter_length': 0, 'max_reporter_length': 0,
            'organism_counts': Counter(), 'full_organism_counts': Counter(),
            'virus_name_counts': Counter(), 'virus_genus_counts': Counter(),
            'virus_family_counts': Counter(), 'genome_type_counts': Counter(),
            'strains_per_organism': {},
            'matched_genomes': Counter(),
            'matched_genomes_by_id': Counter(),
            'matched_genomes_by_name': Counter(),
            'total_distinct_genomes': 0
        }
    
    print(f"Found {len(valid_designs)} valid designs with complete probe information")
    
    # Basic counts
    stats['total_designs'] = len(valid_designs)
    
    # Extract data for easier processing
    capture_gc = [d.get('capture_gc', 0) for d in valid_designs]
    reporter_gc = [d.get('reporter_gc', 0) for d in valid_designs]
    capture_tm = [d.get('capture_tm', 0) for d in valid_designs]
    reporter_tm = [d.get('reporter_tm', 0) for d in valid_designs]
    
    # GC content statistics
    stats['avg_capture_gc'] = np.mean(capture_gc)
    stats['std_capture_gc'] = np.std(capture_gc)
    stats['min_capture_gc'] = min(capture_gc)
    stats['max_capture_gc'] = max(capture_gc)
    
    stats['avg_reporter_gc'] = np.mean(reporter_gc)
    stats['std_reporter_gc'] = np.std(reporter_gc)
    stats['min_reporter_gc'] = min(reporter_gc)
    stats['max_reporter_gc'] = max(reporter_gc)
    
    # Melting temperature statistics
    stats['avg_capture_tm'] = np.mean(capture_tm)
    stats['std_capture_tm'] = np.std(capture_tm)
    stats['min_capture_tm'] = min(capture_tm)
    stats['max_capture_tm'] = max(capture_tm)
    
    stats['avg_reporter_tm'] = np.mean(reporter_tm)
    stats['std_reporter_tm'] = np.std(reporter_tm)
    stats['min_reporter_tm'] = min(reporter_tm)
    stats['max_reporter_tm'] = max(reporter_tm)
    
    # Sequence length statistics - PROPERLY CALCULATED
    target_lengths = [len(d.get('target_sequence', '')) for d in valid_designs]
    capture_lengths = [len(d.get('capture_probe', '')) for d in valid_designs]
    reporter_lengths = [len(d.get('reporter_probe', '')) for d in valid_designs]
    
    stats['avg_target_length'] = np.mean(target_lengths)
    stats['std_target_length'] = np.std(target_lengths)
    stats['min_target_length'] = min(target_lengths)
    stats['max_target_length'] = max(target_lengths)
    
    stats['avg_capture_length'] = np.mean(capture_lengths)
    stats['std_capture_length'] = np.std(capture_lengths)
    stats['min_capture_length'] = min(capture_lengths)
    stats['max_capture_length'] = max(capture_lengths)
    
    stats['avg_reporter_length'] = np.mean(reporter_lengths)
    stats['std_reporter_length'] = np.std(reporter_lengths)
    stats['min_reporter_length'] = min(reporter_lengths)
    stats['max_reporter_length'] = max(reporter_lengths)
    
    print(f"Sequence length summary:")
    print(f"  Target: {stats['avg_target_length']:.2f} ± {stats['std_target_length']:.2f} bases (range: {stats['min_target_length']}-{stats['max_target_length']})")
    print(f"  Capture probe: {stats['avg_capture_length']:.2f} ± {stats['std_capture_length']:.2f} bases (range: {stats['min_capture_length']}-{stats['max_capture_length']})")
    print(f"  Reporter probe: {stats['avg_reporter_length']:.2f} ± {stats['std_reporter_length']:.2f} bases (range: {stats['min_reporter_length']}-{stats['max_reporter_length']})")
    
    # More detailed virus statistics that properly deal with strains
    # Group by main organism name (without strain information)
    organisms = [d.get('organism', 'Unknown') for d in valid_designs]
    full_organism_infos = [d.get('full_organism_info', 'Unknown') for d in valid_designs]
    virus_names = [d.get('virus_name', 'Unknown') for d in valid_designs]
    virus_genera = [d.get('virus_genus', 'Unknown') for d in valid_designs]
    virus_families = [d.get('virus_family', 'Unknown') for d in valid_designs]
    genome_types = [d.get('genome_type', 'Unknown') for d in valid_designs]
    
    # Store the main organism stats (base organism without strain)
    stats['organism_counts'] = Counter(organisms)
    # Store detailed strain information
    stats['full_organism_counts'] = Counter(full_organism_infos)
    stats['virus_name_counts'] = Counter(virus_names)
    stats['virus_genus_counts'] = Counter(virus_genera)
    stats['virus_family_counts'] = Counter(virus_families)
    stats['genome_type_counts'] = Counter(genome_types)
    
    # Count the strains per organism
    organism_to_strains = {}
    for d in valid_designs:
        org = d.get('organism', 'Unknown')
        full_info = d.get('full_organism_info', 'Unknown')
        if org not in organism_to_strains:
            organism_to_strains[org] = set()
        organism_to_strains[org].add(full_info)
    
    # Count the number of strains per organism
    stats['strains_per_organism'] = {org: len(strains) for org, strains in organism_to_strains.items()}
    
    # Match designs to virus genomes in the aligned_mammal_virus_genomes folder
    if genbank_id_map is not None and virus_name_map is not None:
        matched_genomes = set()  # Set of file paths to matched genome files
        matched_by_id = Counter()  # Count matches by GenBank ID
        matched_by_name = Counter()  # Count matches by virus name
        
        # Track designs that match to each genome file
        genome_file_to_designs = {}
        
        for design in valid_designs:
            genbank_id = design.get('genbank_id', '')
            virus_name = design.get('virus_name', '')
            design_matched = False
            
            # Try to match by GenBank ID
            if genbank_id and genbank_id in genbank_id_map:
                genome_file = genbank_id_map[genbank_id]
                matched_genomes.add(genome_file)
                matched_by_id[genome_file] += 1
                design_matched = True
                
                # Add this design to the set of designs for this genome file
                if genome_file not in genome_file_to_designs:
                    genome_file_to_designs[genome_file] = set()
                genome_file_to_designs[genome_file].add(genbank_id)
            
            # If not matched by ID, try to match by virus name
            if not design_matched and virus_name and virus_name in virus_name_map:
                genome_file = virus_name_map[virus_name]
                matched_genomes.add(genome_file)
                matched_by_name[genome_file] += 1
                design_matched = True
                
                # Add this design to the set of designs for this genome file
                if genome_file not in genome_file_to_designs:
                    genome_file_to_designs[genome_file] = set()
                genome_file_to_designs[genome_file].add(virus_name)
        
        # Store the results in stats
        stats['matched_genomes'] = Counter({os.path.basename(file): 1 for file in matched_genomes})
        stats['matched_genomes_by_id'] = Counter({os.path.basename(file): count for file, count in matched_by_id.items()})
        stats['matched_genomes_by_name'] = Counter({os.path.basename(file): count for file, count in matched_by_name.items()})
        stats['total_distinct_genomes'] = len(matched_genomes)
        
        # Print summary of genome matching
        print(f"Genome matching summary:")
        print(f"  Total distinct genomes matched: {len(matched_genomes)}")
        print(f"  Matches by GenBank ID: {sum(matched_by_id.values())}")
        print(f"  Matches by virus name: {sum(matched_by_name.values())}")
        
        # Debug: Print out the matches for verification
        print("\nMatched genome files:")
        for genome_file, design_ids in genome_file_to_designs.items():
            print(f"  {os.path.basename(genome_file)}: {len(design_ids)} designs")
    
    return stats

def format_tables(stats):
    """
    Format statistics into pandas DataFrames for tabular display.
    
    Args:
        stats (dict): Dictionary containing calculated statistics
        
    Returns:
        dict: Dictionary containing formatted tables
    """
    tables = {}
    
    # Basic statistics table with length information
    basic_stats = {
        'Statistic': [
            'Total Valid Designs',
            'Target Length (avg)',
            'Target Length (min-max)',
            'Capture Probe Length (avg)',
            'Capture Probe Length (min-max)',
            'Reporter Probe Length (avg)',
            'Reporter Probe Length (min-max)'
        ],
        'Value': [
            stats['total_designs'],
            f"{stats['avg_target_length']:.2f} ± {stats['std_target_length']:.2f}",
            f"{stats['min_target_length']} - {stats['max_target_length']}",
            f"{stats['avg_capture_length']:.2f} ± {stats['std_capture_length']:.2f}",
            f"{stats['min_capture_length']} - {stats['max_capture_length']}",
            f"{stats['avg_reporter_length']:.2f} ± {stats['std_reporter_length']:.2f}",
            f"{stats['min_reporter_length']} - {stats['max_reporter_length']}"
        ]
    }
    tables['basic_stats'] = pd.DataFrame(basic_stats)
    
    # GC content statistics table
    gc_stats = {
        'Probe Type': ['Capture Probe', 'Reporter Probe'],
        'Avg GC (%)': [
            f"{stats['avg_capture_gc']:.2f}",
            f"{stats['avg_reporter_gc']:.2f}"
        ],
        'Std Dev': [
            f"{stats['std_capture_gc']:.2f}",
            f"{stats['std_reporter_gc']:.2f}"
        ],
        'Min (%)': [
            f"{stats['min_capture_gc']:.2f}",
            f"{stats['min_reporter_gc']:.2f}"
        ],
        'Max (%)': [
            f"{stats['max_capture_gc']:.2f}",
            f"{stats['max_reporter_gc']:.2f}"
        ]
    }
    tables['gc_stats'] = pd.DataFrame(gc_stats)
    
    # Melting temperature statistics table
    tm_stats = {
        'Probe Type': ['Capture Probe', 'Reporter Probe'],
        'Avg Tm (°C)': [
            f"{stats['avg_capture_tm']:.2f}",
            f"{stats['avg_reporter_tm']:.2f}"
        ],
        'Std Dev': [
            f"{stats['std_capture_tm']:.2f}",
            f"{stats['std_reporter_tm']:.2f}"
        ],
        'Min (°C)': [
            f"{stats['min_capture_tm']:.2f}",
            f"{stats['min_reporter_tm']:.2f}"
        ],
        'Max (°C)': [
            f"{stats['max_capture_tm']:.2f}",
            f"{stats['max_reporter_tm']:.2f}"
        ]
    }
    tables['tm_stats'] = pd.DataFrame(tm_stats)
    
    # Organism distribution (without strain info)
    organism_data = []
    for org, count in stats['organism_counts'].most_common():
        strain_count = stats['strains_per_organism'].get(org, 0)
        organism_data.append({
            'Organism': org,
            'Count': count,
            'Percent': count / stats['total_designs'] * 100 if stats['total_designs'] > 0 else 0,
            'Strain Count': strain_count
        })
    tables['organism_stats'] = pd.DataFrame(organism_data)
    
    # Matched genome statistics
    if 'matched_genomes' in stats and stats['matched_genomes']:
        genome_data = []
        for genome_file, count in stats['matched_genomes'].most_common():
            by_id = stats['matched_genomes_by_id'].get(genome_file, 0)
            by_name = stats['matched_genomes_by_name'].get(genome_file, 0)
            genome_data.append({
                'Genome File': genome_file,
                'Matched Designs': by_id + by_name,
                'Matched by GenBank ID': by_id,
                'Matched by Virus Name': by_name
            })
        tables['genome_stats'] = pd.DataFrame(genome_data)
    
    # Virus family distribution table
    virus_family_data = [
        {'Family': family, 'Count': count, 'Percentage': (count / stats['total_designs']) * 100}
        for family, count in stats['virus_family_counts'].most_common(30)  # Top 30 families
    ]
    tables['virus_family'] = pd.DataFrame(virus_family_data)
    
    # Virus genus distribution table
    virus_genus_data = [
        {'Genus': genus, 'Count': count, 'Percentage': (count / stats['total_designs']) * 100}
        for genus, count in stats['virus_genus_counts'].most_common(10)  # Top 10 genera
    ]
    tables['virus_genus'] = pd.DataFrame(virus_genus_data)
    
    # Genome type distribution table
    genome_type_data = [
        {'Genome Type': genome_type, 'Count': count, 'Percentage': (count / stats['total_designs']) * 100}
        for genome_type, count in stats['genome_type_counts'].most_common(20)  # Top 20 genome types
    ]
    tables['genome_type'] = pd.DataFrame(genome_type_data)
    
    # Top strain data (for reference)
    strain_data = [
        {'Full Organism Info': strain, 'Count': count, 'Percentage': (count / stats['total_designs']) * 100}
        for strain, count in stats['full_organism_counts'].most_common(30)  # Top 30 strains
    ]
    tables['strain_stats'] = pd.DataFrame(strain_data)
    
    return tables

def save_tables_to_file(tables, output_file):
    """
    Save formatted tables to a text file.
    
    Args:
        tables (dict): Dictionary containing formatted tables
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        f.write("# QDB Assay Design Statistics\n\n")
        
        f.write("## Basic Statistics\n")
        f.write(tabulate(tables['basic_stats'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## GC Content Statistics\n")
        f.write(tabulate(tables['gc_stats'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Melting Temperature Statistics\n")
        f.write(tabulate(tables['tm_stats'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Top 30 Organisms (Base Name)\n")
        f.write(tabulate(tables['organism_stats'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Top 30 Virus Strains\n")
        f.write(tabulate(tables['strain_stats'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Top 30 Virus Families\n")
        f.write(tabulate(tables['virus_family'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Top 10 Virus Genera\n")
        f.write(tabulate(tables['virus_genus'], headers='keys', tablefmt='grid'))
        f.write("\n\n")
        
        f.write("## Genome Type Distribution\n")
        f.write(tabulate(tables['genome_type'], headers='keys', tablefmt='grid'))
        
        if 'genome_stats' in tables:
            f.write("\n\n## Matched Genome Statistics\n")
            f.write(tabulate(tables['genome_stats'], headers='keys', tablefmt='grid'))

def save_tables_to_csv(tables, output_prefix):
    """
    Save formatted tables to CSV files.
    
    Args:
        tables (dict): Dictionary containing formatted tables
        output_prefix (str): Prefix for the output file names
    """
    # Save basic statistics
    tables['basic_stats'].to_csv(f"{output_prefix}_basic_stats.csv", index=False)
    
    # Save GC content statistics
    tables['gc_stats'].to_csv(f"{output_prefix}_gc_stats.csv", index=False)
    
    # Save melting temperature statistics
    tables['tm_stats'].to_csv(f"{output_prefix}_tm_stats.csv", index=False)
    
    # Save organism distribution
    tables['organism_stats'].to_csv(f"{output_prefix}_organism_stats.csv", index=False)
    
    # Save strain statistics
    tables['strain_stats'].to_csv(f"{output_prefix}_strain_stats.csv", index=False)
    
    # Save virus family distribution
    tables['virus_family'].to_csv(f"{output_prefix}_virus_family.csv", index=False)
    
    # Save virus genus distribution
    tables['virus_genus'].to_csv(f"{output_prefix}_virus_genus.csv", index=False)
    
    # Save genome type distribution
    tables['genome_type'].to_csv(f"{output_prefix}_genome_type.csv", index=False)
    
    # Save matched genome statistics
    if 'genome_stats' in tables:
        tables['genome_stats'].to_csv(f"{output_prefix}_matched_genome_stats.csv", index=False)
    
    # Save all statistics in a single file
    with open(f"{output_prefix}_all_stats.csv", 'w', newline='') as f:
        f.write("# QDB Assay Design Statistics\n\n")
        
        f.write("## Basic Statistics\n")
        tables['basic_stats'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## GC Content Statistics\n")
        tables['gc_stats'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Melting Temperature Statistics\n")
        tables['tm_stats'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Top 30 Organisms\n")
        tables['organism_stats'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Top 30 Virus Strains\n")
        tables['strain_stats'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Top 30 Virus Families\n")
        tables['virus_family'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Top 10 Virus Genera\n")
        tables['virus_genus'].to_csv(f, index=False)
        f.write("\n")
        
        f.write("## Genome Type Distribution\n")
        tables['genome_type'].to_csv(f, index=False)
        
    print(f"All statistics saved to {output_prefix}_all_stats.csv")

def save_raw_probe_data(designs, output_file):
    """
    Save raw probe data to a CSV file.
    
    Args:
        designs (list): List of dictionaries containing design information
        output_file (str): Path to output file
    """
    # Extract relevant fields for each design
    probe_data = []
    
    for i, design in enumerate(designs):
        # Skip designs without complete probe information
        if not ('target_sequence' in design and 
                'capture_probe' in design and 
                'reporter_probe' in design and
                'capture_gc' in design and
                'reporter_gc' in design):
            continue
            
        # Calculate lengths - ensure we're only counting the actual sequence length
        target_length = len(design.get('target_sequence', ''))
        capture_length = len(design.get('capture_probe', ''))
        reporter_length = len(design.get('reporter_probe', ''))
        
        # Extract data
        entry = {
            'design_id': i + 1,
            'genbank_id': design.get('genbank_id', ''),
            'organism': design.get('organism', ''),
            'full_organism_info': design.get('full_organism_info', ''),
            'virus_name': design.get('virus_name', ''),
            'virus_genus': design.get('virus_genus', ''),
            'virus_family': design.get('virus_family', ''),
            'genome_type': design.get('genome_type', ''),
            'target_sequence': design.get('target_sequence', ''),
            'target_length': target_length,
            'capture_probe': design.get('capture_probe', ''),
            'capture_length': capture_length,
            'capture_gc': design.get('capture_gc', 0),
            'capture_tm': design.get('capture_tm', 0),
            'reporter_probe': design.get('reporter_probe', ''),
            'reporter_length': reporter_length,
            'reporter_gc': design.get('reporter_gc', 0),
            'reporter_tm': design.get('reporter_tm', 0)
        }
        
        probe_data.append(entry)
    
    if not probe_data:
        print("ERROR: No valid probe data to save")
        return
        
    # Convert to DataFrame and save to CSV
    df = pd.DataFrame(probe_data)
    df.to_csv(output_file, index=False)
    
    print(f"Raw probe data saved to {output_file} - total of {len(probe_data)} valid probe designs")

def main():
    """Main function to run the analysis."""
    input_file = "qdb_assays_output_all_genomes.txt"
    output_file = "qdb_design_statistics.txt"
    output_csv_prefix = "qdb_design_statistics"
    output_raw_data = "qdb_design_raw_data.csv"
    aligned_genomes_dir = "aligned_mammal_virus_genomes"
    
    print(f"Parsing design file: {input_file}")
    designs = parse_qdb_design_file(input_file)
    print(f"Found {len(designs)} total designs in file")
    
    # Extract virus genome information from FASTA files
    print(f"Extracting virus genome information from {aligned_genomes_dir}...")
    genbank_id_map, virus_name_map = extract_virus_genomes_from_fasta_directory(aligned_genomes_dir)
    
    print("Calculating statistics for designs...")
    stats = calculate_statistics(designs, genbank_id_map, virus_name_map)
    
    print("Formatting tables...")
    tables = format_tables(stats)
    
    print(f"Saving tables to: {output_file}")
    save_tables_to_file(tables, output_file)
    
    print(f"Saving tables to CSV files with prefix: {output_csv_prefix}")
    save_tables_to_csv(tables, output_csv_prefix)
    
    print(f"Saving raw probe data to: {output_raw_data}")
    save_raw_probe_data(designs, output_raw_data)
    
    print("\nBasic Statistics:")
    print(tabulate(tables['basic_stats'], headers='keys', tablefmt='grid'))
    
    print("\nGC Content Statistics:")
    print(tabulate(tables['gc_stats'], headers='keys', tablefmt='grid'))
    
    print("\nMelting Temperature Statistics:")
    print(tabulate(tables['tm_stats'], headers='keys', tablefmt='grid'))
    
    print("\nSequence Length Summary:")
    print(f"Target sequence: {stats['avg_target_length']:.2f} ± {stats['std_target_length']:.2f} bases")
    print(f"Capture probe: {stats['avg_capture_length']:.2f} ± {stats['std_capture_length']:.2f} bases")
    print(f"Reporter probe: {stats['avg_reporter_length']:.2f} ± {stats['std_reporter_length']:.2f} bases")
    
    print("\nAnalysis complete! Full results saved to:")
    print(f"- Text format: {output_file}")
    print(f"- CSV statistics: {output_csv_prefix}_all_stats.csv")
    print(f"- Raw probe data: {output_raw_data}")

if __name__ == "__main__":
    main()
