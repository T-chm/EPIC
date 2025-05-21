#!/usr/bin/env python3
"""
check_probe_specificity_simplified.py
Simplified version that screens QDB probes and PCR primers against a local virus database
by reading directly from combined design files and outputting specificity results.

Usage:
    python check_probe_specificity_simplified.py
"""

# Import necessary libraries
import subprocess
import tempfile
import textwrap
import os
import re
from pathlib import Path
import time
import math
from Bio import SeqIO
from io import StringIO

# Define constants and configuration
IDENTITY_THRESHOLD = 80.0  # Minimum percent identity to consider as cross-reactive
COV_THRESHOLD = 80.0      # Minimum query coverage to consider as cross-reactive
BLAST_TIMEOUT = 300     # Seconds before we consider a BLAST job to be hanging
MIN_ALIGN_LENGTH = 12    # Minimum alignment length for potential binding site
MAX_3PRIME_MISMATCH = 3  # Maximum allowed mismatches at 3' end
MAX_TOTAL_MISMATCH = 9   # Maximum total mismatches allowed
DELTAG_THRESHOLD = -8.0  # Delta G threshold for binding (kcal/mol)

# BLAST parameters optimized for each sequence type
PARAMS = {
    "primer": dict(task="blastn-short", word_size=7,  evalue="1000"),
    "probe": dict(task="blastn", word_size=11, evalue="1e-5")
}

def calculate_delta_g(sequence1, sequence2):
    """Calculate the delta G (Gibbs free energy) for hybridization between two sequences.
    Uses nearest-neighbor model with empirical parameters for DNA/DNA hybridization."""
    # Nearest-neighbor thermodynamic parameters for DNA/DNA duplexes (SantaLucia, 1998)
    # Values are in kcal/mol
    nn_parameters = {
        'AA': -1.0, 'TT': -1.0,
        'AT': -0.9, 'TA': -0.9,
        'CA': -1.5, 'GT': -1.5,
        'CT': -1.3, 'GA': -1.3,
        'CG': -2.1, 'GC': -2.1,
        'GG': -1.8, 'CC': -1.8
    }
    
    # Ensure sequences are of equal length by truncating the longer one if needed
    min_len = min(len(sequence1), len(sequence2))
    seq1 = sequence1[:min_len].upper()
    seq2 = sequence2[:min_len].upper()
    
    # Reverse complement sequence2 for proper alignment
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    seq2_rc = ''.join(complement.get(base, 'N') for base in reversed(seq2))
    
    # Calculate delta G from nearest-neighbor model
    delta_g = 0.0
    for i in range(len(seq1) - 1):
        pair1 = seq1[i:i+2]
        pair2 = seq2_rc[i:i+2]
        
        # Skip calculation for any pair containing non-standard bases
        if 'N' in pair1 or 'N' in pair2:
            continue
            
        # Calculate mismatch penalty
        if seq1[i] != seq2_rc[i] or seq1[i+1] != seq2_rc[i+1]:
            # Mismatches penalize binding
            delta_g += 0.5  # Penalty for mismatches
        else:
            # Proper base pairing contributes to binding energy
            delta_g += nn_parameters.get(pair1, 0)
    
    # Add initiation penalty
    delta_g += 0.2
    
    # Temperature adjustment (assuming 37°C)
    delta_g += 0.0027 * (310.15 - 273.15) * min_len
    
    return delta_g

def get_sequence_from_db(blast_db, subject_id, start, end):
    """Extract a sequence fragment from the BLAST database."""
    if start > end:
        start, end = end, start  # Ensure start < end
        
    cmd = [
        "blastdbcmd",
        "-db", blast_db,
        "-entry", subject_id,
        "-range", f"{start}-{end}"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Parse the FASTA output to extract just the sequence
        fasta = result.stdout.strip()
        if fasta:
            seq_record = SeqIO.read(StringIO(fasta), "fasta")
            return str(seq_record.seq)
        return ""
    except Exception as e:
        print(f"Error extracting sequence: {e}")
        return ""

def count_3prime_mismatches(query_seq, subject_seq, length=5):
    """Count mismatches at the 3' end of the query sequence."""
    # Ensure sequences are of equal length
    min_len = min(len(query_seq), len(subject_seq))
    query = query_seq[:min_len]
    subject = subject_seq[:min_len]
    
    # Get the 3' end (last few bases) of the query
    end_length = min(length, min_len)
    query_end = query[-end_length:]
    subject_end = subject[-end_length:]
    
    # Count mismatches
    mismatches = sum(q != s for q, s in zip(query_end, subject_end))
    return mismatches

def run_blast(fasta_file, output_file, db_path, seq_type):
    """Run BLAST search for the sequences in the provided FASTA file."""
    p = PARAMS[seq_type]
    
    # Updated BLAST command to output sequence data
    cmd = [
        "blastn", 
        "-task", p["task"], 
        "-query", fasta_file, 
        "-db", db_path,
        "-word_size", str(p["word_size"]),
        "-dust", "no", 
        "-soft_masking", "false",
        "-penalty", "-2", 
        "-reward", "1",
        "-gapopen", "5", 
        "-gapextend", "2",
        "-strand", "both", 
        "-evalue", p["evalue"],
        "-max_target_seqs", "10",  # Increased to get more potential hits
        # Using custom output format to get query and subject sequences
        "-outfmt", "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore qseq sseq",
        "-num_threads", "2", 
        "-out", output_file
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, 
                      stderr=subprocess.DEVNULL, timeout=BLAST_TIMEOUT)
        return True
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
        print(f"Warning: BLAST error or timeout for {fasta_file}")
        return False

def identify_cross_reactive(blast_result_file, blast_db, seq_type):
    """Identify cross-reactive sequences using the new primer binding criteria."""
    if not os.path.exists(blast_result_file) or os.path.getsize(blast_result_file) == 0:
        return set(), {}
    
    bad_ids = set()
    blast_details = {}
    primers = {}
    
    # Parse the BLAST results
    with open(blast_result_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 13:  # Need at least 13 fields with our custom outfmt
                continue
                
            qseqid = parts[0]     # Query ID
            sseqid = parts[1]     # Subject ID
            pident = float(parts[2])  # Percent identity
            length = int(parts[3])    # Alignment length
            mismatches = int(parts[4])  # Number of mismatches
            qstart = int(parts[5])    # Query start position
            qend = int(parts[6])      # Query end position
            sstart = int(parts[7])    # Subject start position
            send = int(parts[8])      # Subject end position
            evalue = float(parts[9])  # E-value
            bitscore = float(parts[10])  # Bit score
            qseq = parts[11]     # Query sequence
            sseq = parts[12]     # Subject sequence
            
            # Skip self-hits
            if qseqid == sseqid:
                continue
                
            try:
                # Store the primer sequence if we haven't seen it yet
                if qseqid not in primers and qseq:
                    primers[qseqid] = qseq
                
                # Initialize blast details for this query
                if qseqid not in blast_details:
                    blast_details[qseqid] = []
                    
                # Calculate delta G only for primer-like sequences
                delta_g = None
                three_prime_mismatches = None
                is_binding_site = False
                
                # Apply the four criteria for potential binding sites
                if seq_type == "primer":
                    # 1) Aligned genome sequence > 12 bp and e-value criteria
                    if length >= MIN_ALIGN_LENGTH:
                        # 2) Calculate delta G
                        delta_g = calculate_delta_g(qseq, sseq)
                        
                        # 3) Check 3' end mismatches
                        three_prime_mismatches = count_3prime_mismatches(qseq, sseq)
                        
                        # 4) Check total mismatches
                        is_binding_site = (
                            length >= MIN_ALIGN_LENGTH and
                            (evalue > 1000 or evalue < 1e-10) and  # Either very high or very low e-value
                            delta_g <= DELTAG_THRESHOLD and
                            three_prime_mismatches < MAX_3PRIME_MISMATCH and
                            mismatches < MAX_TOTAL_MISMATCH
                        )
                        
                        if is_binding_site:
                            bad_ids.add(qseqid)
                else:
                    # For probes, use the original identity/coverage criteria
                    coverage = 100 * length / (qend - qstart + 1)
                    is_binding_site = pident >= IDENTITY_THRESHOLD and coverage >= COV_THRESHOLD
                    if is_binding_site:
                        bad_ids.add(qseqid)
                
                # Store detailed hit information
                if len(blast_details[qseqid]) < 3 or is_binding_site:  # Always include binding sites
                    hit_info = {
                        'subject_id': sseqid,
                        'percent_identity': pident,
                        'alignment_length': length,
                        'mismatches': mismatches,
                        'e_value': evalue,
                        'bit_score': bitscore,
                        'query_seq': qseq,
                        'subject_seq': sseq,
                        'is_cross_reactive': is_binding_site
                    }
                    
                    # Add primer-specific details if available
                    if delta_g is not None:
                        hit_info['delta_g'] = delta_g
                    if three_prime_mismatches is not None:
                        hit_info['three_prime_mismatches'] = three_prime_mismatches
                    
                    # Replace less significant hits with this one if it's a binding site
                    if is_binding_site and len(blast_details[qseqid]) >= 3:
                        # Find a non-binding site to replace
                        for i, existing_hit in enumerate(blast_details[qseqid]):
                            if not existing_hit['is_cross_reactive']:
                                blast_details[qseqid][i] = hit_info
                                break
                        # If all are binding sites, just append (we'll sort later)
                        if all(hit['is_cross_reactive'] for hit in blast_details[qseqid]):
                            blast_details[qseqid].append(hit_info)
                    else:
                        blast_details[qseqid].append(hit_info)
                    
            except (ValueError, IndexError, ZeroDivisionError) as e:
                print(f"Error processing BLAST hit: {e}")
                continue
    
    # Sort hits for each query by cross-reactivity and e-value
    for qseqid in blast_details:
        blast_details[qseqid].sort(key=lambda x: (not x['is_cross_reactive'], x['e_value']))
        # Limit to top 3 hits after sorting
        blast_details[qseqid] = blast_details[qseqid][:3]
    
    return bad_ids, blast_details

def parse_combined_designs(file_path):
    """Parse the combined QDB probe and PCR primer designs from a single file."""
    designs = []
    current_design = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith("Target Sequence:"):
            # Start of a new design
            if "genbank_id" in current_design and current_design["genbank_id"]:
                designs.append(current_design)
                current_design = {}
            
            current_design["target_seq"] = line.replace("Target Sequence:", "").strip()
        
        elif line.startswith("Capture Probe:"):
            # Extract sequence between /5AmMC6/ and space
            match = re.search(r'/5AmMC6/([ATCG]+)', line)
            if match:
                current_design["capture_probe"] = match.group(1)
        
        elif line.startswith("Reporter Probe:"):
            # Extract sequence between beginning and /3Cy5Sp/
            match = re.search(r'^Reporter Probe:\s+([ATCG]+)/3Cy5Sp/', line)
            if match:
                current_design["reporter_probe"] = match.group(1)
        
        elif line.startswith("GenBank ID:"):
            current_design["genbank_id"] = line.replace("GenBank ID:", "").strip()
        
        elif line.startswith("Metadata:"):
            current_design["metadata"] = line.replace("Metadata:", "").strip()
        
        elif line.startswith("Forward Primer:"):
            # Extract sequence before first space or parenthesis
            match = re.search(r'^Forward Primer:\s+([ATCG]+)', line)
            if match:
                current_design["forward_primer"] = match.group(1)
        
        elif line.startswith("Reverse Primer:"):
            # Extract sequence before first space or parenthesis
            match = re.search(r'^Reverse Primer:\s+([ATCG]+)', line)
            if match:
                current_design["reverse_primer"] = match.group(1)
        
        i += 1
    
    # Add the last design if it exists
    if "genbank_id" in current_design and current_design["genbank_id"]:
        designs.append(current_design)
    
    return designs

def create_sequence_dict(combined_designs):
    """Extract all sequences from the combined designs into a dictionary."""
    sequences = {}
    
    for i, design in enumerate(combined_designs):
        genbank_id = design.get("genbank_id", "")
        if not genbank_id:
            continue
            
        # Extract QDB probes
        if "capture_probe" in design:
            seq_id = f"{genbank_id}_CAP_{i}"
            sequences[seq_id] = design["capture_probe"]
            
        if "reporter_probe" in design:
            seq_id = f"{genbank_id}_REP_{i}"
            sequences[seq_id] = design["reporter_probe"]
            
        # Extract PCR primers
        if "forward_primer" in design:
            seq_id = f"{genbank_id}_F_{i}"
            sequences[seq_id] = design["forward_primer"]
            
        if "reverse_primer" in design:
            seq_id = f"{genbank_id}_R_{i}"
            sequences[seq_id] = design["reverse_primer"]
    
    return sequences

def write_sequences_to_fasta(sequences, output_file):
    """Write sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n{textwrap.fill(sequence, 80)}\n")

def check_sequence_specificity(sequences, blast_db, temp_dir, seq_type):
    """Check specificity of sequences by running BLAST."""
    # Create temporary FASTA file
    fasta_file = os.path.join(temp_dir, f"{seq_type}_sequences.fa")
    write_sequences_to_fasta(sequences, fasta_file)
    
    # Run BLAST
    blast_out = os.path.join(temp_dir, f"{seq_type}_blast_results.tsv")
    success = run_blast(fasta_file, blast_out, blast_db, seq_type)
    
    # Process BLAST results
    bad_ids = set()
    blast_details = {}
    if success:
        bad_ids, blast_details = identify_cross_reactive(blast_out, blast_db, seq_type)
    
    return bad_ids, blast_details

def write_results_with_specificity(combined_designs, good_ids, blast_details, output_path):
    """Write all designs with specificity annotations and detailed BLAST results to output file."""
    # Track count of specific designs
    fully_specific_count = 0
    total_designs = 0
    
    # Write all designs with specificity annotations
    with open(output_path, 'w') as f:
        # Write header with explanation
        f.write("# BLAST Specificity Check Results\n")
        f.write("# All designs are shown with specificity check results for each component\n")
        f.write("# [SPECIFIC] - Component passed specificity check (no significant cross-reactivity)\n")
        f.write("# [CROSS-REACTIVE] - Component failed specificity check (has significant cross-reactivity)\n\n")
        
        # Write primer binding criteria
        f.write("# Primer Binding Site Criteria:\n")
        f.write(f"#  1) Aligned sequence length ≥ {MIN_ALIGN_LENGTH} bp\n")
        f.write(f"#  2) Delta G value ≤ {DELTAG_THRESHOLD} kcal/mol\n")
        f.write(f"#  3) 3' end mismatches < {MAX_3PRIME_MISMATCH} bp\n")
        f.write(f"#  4) Total mismatches < {MAX_TOTAL_MISMATCH} bp\n\n")
        
        # Write probe criteria
        f.write("# Probe Cross-Reactivity Criteria:\n")
        f.write(f"#  - Identity threshold: {IDENTITY_THRESHOLD}%\n")
        f.write(f"#  - Coverage threshold: {COV_THRESHOLD}%\n\n")
        
        for i, design in enumerate(combined_designs):
            genbank_id = design.get("genbank_id", "")
            if not genbank_id:
                continue
                
            total_designs += 1
            
            # Track specificity for each component
            qdb_cap_specific = True
            qdb_rep_specific = True
            pcr_fwd_specific = True
            pcr_rev_specific = True
            
            # Sequence IDs for this design
            cap_id = f"{genbank_id}_CAP_{i}"
            rep_id = f"{genbank_id}_REP_{i}"
            fwd_id = f"{genbank_id}_F_{i}"
            rev_id = f"{genbank_id}_R_{i}"
            
            # Check QDB probe specificity
            if "capture_probe" in design and "reporter_probe" in design:
                qdb_cap_specific = cap_id in good_ids
                qdb_rep_specific = rep_id in good_ids
            
            # Check PCR primer specificity
            if "forward_primer" in design and "reverse_primer" in design:
                pcr_fwd_specific = fwd_id in good_ids
                pcr_rev_specific = rev_id in good_ids
            
            # Check if all components are specific
            all_specific = (qdb_cap_specific and qdb_rep_specific and 
                          pcr_fwd_specific and pcr_rev_specific)
            
            if all_specific:
                fully_specific_count += 1
            
            # Write target sequence if available
            if "target_seq" in design:
                f.write(f"Target Sequence: {design['target_seq']}\n")
            
            # Write QDB probe info with specificity status and detailed BLAST results
            if "capture_probe" in design and "reporter_probe" in design:
                cap_status = "[SPECIFIC]" if qdb_cap_specific else "[CROSS-REACTIVE]"
                rep_status = "[SPECIFIC]" if qdb_rep_specific else "[CROSS-REACTIVE]"
                
                f.write(f"Capture Probe: /5AmMC6/{design['capture_probe']} {cap_status}\n")
                # Write detailed BLAST results for capture probe if available
                if cap_id in blast_details and blast_details[cap_id]:
                    f.write("  BLAST Hits:\n")
                    for hit in blast_details[cap_id]:
                        status = "[CROSS-REACTIVE]" if hit['is_cross_reactive'] else ""
                        f.write(f"  - Hit: {hit['subject_id']} ")
                        f.write(f"Identity: {hit['percent_identity']:.1f}% ")
                        
                        # Write coverage if available (for probes)
                        if 'coverage' in hit:
                            f.write(f"Coverage: {hit['coverage']:.1f}% ")
                            
                        f.write(f"Mismatches: {hit['mismatches']} ")
                        f.write(f"E-value: {hit['e_value']:.2e} {status}\n")
                        
                        # Show sequences for cross-reactive hits
                        if hit['is_cross_reactive'] and 'query_seq' in hit and 'subject_seq' in hit:
                            f.write(f"    Query: {hit['query_seq']}\n")
                            f.write(f"    Subj:  {hit['subject_seq']}\n")
                
                f.write(f"Reporter Probe: {design['reporter_probe']}/3Cy5Sp/ {rep_status}\n")
                # Write detailed BLAST results for reporter probe if available
                if rep_id in blast_details and blast_details[rep_id]:
                    f.write("  BLAST Hits:\n")
                    for hit in blast_details[rep_id]:
                        status = "[CROSS-REACTIVE]" if hit['is_cross_reactive'] else ""
                        f.write(f"  - Hit: {hit['subject_id']} ")
                        f.write(f"Identity: {hit['percent_identity']:.1f}% ")
                        
                        # Write coverage if available (for probes)
                        if 'coverage' in hit:
                            f.write(f"Coverage: {hit['coverage']:.1f}% ")
                            
                        f.write(f"Mismatches: {hit['mismatches']} ")
                        f.write(f"E-value: {hit['e_value']:.2e} {status}\n")
                        
                        # Show sequences for cross-reactive hits
                        if hit['is_cross_reactive'] and 'query_seq' in hit and 'subject_seq' in hit:
                            f.write(f"    Query: {hit['query_seq']}\n")
                            f.write(f"    Subj:  {hit['subject_seq']}\n")
            
            # Write GenBank ID and metadata
            f.write(f"GenBank ID: {genbank_id}\n")
            if "metadata" in design:
                f.write(f"Metadata: {design['metadata']}\n")
            
            # Write PCR primer info with specificity status and detailed BLAST results
            if "forward_primer" in design and "reverse_primer" in design:
                fwd_status = "[SPECIFIC]" if pcr_fwd_specific else "[CROSS-REACTIVE]"
                rev_status = "[SPECIFIC]" if pcr_rev_specific else "[CROSS-REACTIVE]"
                
                f.write(f"Forward Primer: {design['forward_primer']} {fwd_status}\n")
                # Write detailed BLAST results for forward primer if available
                if fwd_id in blast_details and blast_details[fwd_id]:
                    f.write("  BLAST Hits:\n")
                    for hit in blast_details[fwd_id]:
                        status = "[CROSS-REACTIVE]" if hit['is_cross_reactive'] else ""
                        f.write(f"  - Hit: {hit['subject_id']} ")
                        f.write(f"Identity: {hit['percent_identity']:.1f}% ")
                        f.write(f"Len: {hit['alignment_length']} bp ")
                        f.write(f"Mismatches: {hit['mismatches']} ")
                        
                        # Write primer-specific criteria
                        if 'delta_g' in hit:
                            f.write(f"Delta G: {hit['delta_g']:.1f} kcal/mol ")
                        if 'three_prime_mismatches' in hit:
                            f.write(f"3' Mismatches: {hit['three_prime_mismatches']} ")
                            
                        f.write(f"E-value: {hit['e_value']:.2e} {status}\n")
                        
                        # Show sequences for cross-reactive hits
                        if hit['is_cross_reactive'] and 'query_seq' in hit and 'subject_seq' in hit:
                            f.write(f"    Query: {hit['query_seq']}\n")
                            f.write(f"    Subj:  {hit['subject_seq']}\n")
                
                f.write(f"Reverse Primer: {design['reverse_primer']} {rev_status}\n")
                # Write detailed BLAST results for reverse primer if available
                if rev_id in blast_details and blast_details[rev_id]:
                    f.write("  BLAST Hits:\n")
                    for hit in blast_details[rev_id]:
                        status = "[CROSS-REACTIVE]" if hit['is_cross_reactive'] else ""
                        f.write(f"  - Hit: {hit['subject_id']} ")
                        f.write(f"Identity: {hit['percent_identity']:.1f}% ")
                        f.write(f"Len: {hit['alignment_length']} bp ")
                        f.write(f"Mismatches: {hit['mismatches']} ")
                        
                        # Write primer-specific criteria
                        if 'delta_g' in hit:
                            f.write(f"Delta G: {hit['delta_g']:.1f} kcal/mol ")
                        if 'three_prime_mismatches' in hit:
                            f.write(f"3' Mismatches: {hit['three_prime_mismatches']} ")
                            
                        f.write(f"E-value: {hit['e_value']:.2e} {status}\n")
                        
                        # Show sequences for cross-reactive hits
                        if hit['is_cross_reactive'] and 'query_seq' in hit and 'subject_seq' in hit:
                            f.write(f"    Query: {hit['query_seq']}\n")
                            f.write(f"    Subj:  {hit['subject_seq']}\n")
            
            # Write amplicon if available
            if "amplicon_sequence" in design:
                f.write(f"Amplicon Sequence: {design['amplicon_sequence']}\n")
            
            # Add overall specificity status
            overall_status = "[ALL SPECIFIC]" if all_specific else "[HAS CROSS-REACTIVITY]"
            f.write(f"Overall Specificity: {overall_status}\n\n")
    
    return total_designs, fully_specific_count

def main():
    t0 = time.time()
    
    # Hardcoded file paths
    combined_file = "new_pcr_primers_output_rep_copy.txt"
    blast_db = "local_virus_blast_db/ref_viruses_rep_genomes"
    
    print(f"Starting specificity check with enhanced primer binding analysis...")
    print(f"Using the following criteria for primer binding:")
    print(f"  1) Aligned sequence length ≥ {MIN_ALIGN_LENGTH} bp")
    print(f"  2) Delta G ≤ {DELTAG_THRESHOLD} kcal/mol")
    print(f"  3) 3' end mismatches < {MAX_3PRIME_MISMATCH} bp")
    print(f"  4) Total mismatches < {MAX_TOTAL_MISMATCH} bp")
    
    # Set up temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # Parse combined designs
        t1 = time.time()
        print(f"\nParsing combined designs from {combined_file}...")
        combined_designs = parse_combined_designs(combined_file)
        print(f"Found {len(combined_designs)} combined design sets in {time.time()-t1:.1f}s")
        
        # Extract all sequences
        all_sequences = create_sequence_dict(combined_designs)
        print(f"Extracted {len(all_sequences)} total sequences to check")
        
        # Separate QDB probes and PCR primers for appropriate BLAST parameters
        qdb_sequences = {seq_id: seq for seq_id, seq in all_sequences.items() 
                        if '_CAP_' in seq_id or '_REP_' in seq_id}
        pcr_sequences = {seq_id: seq for seq_id, seq in all_sequences.items() 
                        if '_F_' in seq_id or '_R_' in seq_id}
        
        # Combined BLAST details dictionary to store all hit information
        all_blast_details = {}
        
        # Process QDB probes
        print(f"\nProcessing {len(qdb_sequences)} QDB probe sequences...")
        qdb_bad_ids, qdb_blast_details = check_sequence_specificity(qdb_sequences, blast_db, temp_dir, "probe")
        all_blast_details.update(qdb_blast_details)
        print(f"Found {len(qdb_bad_ids)} cross-reactive QDB probe sequences")
        
        # Process PCR primers
        print(f"\nProcessing {len(pcr_sequences)} PCR primer sequences using binding site criteria...")
        pcr_bad_ids, pcr_blast_details = check_sequence_specificity(pcr_sequences, blast_db, temp_dir, "primer")
        all_blast_details.update(pcr_blast_details)
        print(f"Found {len(pcr_bad_ids)} PCR primers with potential binding sites")
        
        # Combine cross-reactive IDs
        bad_ids = qdb_bad_ids.union(pcr_bad_ids)
        
        # Determine good IDs (non-cross-reactive sequences)
        good_ids = set(all_sequences.keys()) - bad_ids
        
        # Write results with specificity annotations and detailed BLAST information
        results_out_path = Path(combined_file).with_suffix(".binding_analysis.txt")
        print(f"\nWriting results with detailed binding site analysis to output file...")
        total_designs, fully_specific_count = write_results_with_specificity(
            combined_designs, good_ids, all_blast_details, results_out_path)
        
        # Print summary statistics
        kept = len(good_ids)
        discarded = len(bad_ids)
        total = len(all_sequences)
        pct = 100 * discarded / total if total else 0
        
        print(f"\nAnalysis complete in {time.time()-t0:.1f} seconds")
        print(f"Sequence-level results:")
        print(f"  - Specific: {kept}/{total} sequences ({100-pct:.1f}%)")
        print(f"  - Cross-reactive/binding: {discarded}/{total} sequences ({pct:.1f}%)")
        
        # Show summary of design specificity
        designs_pct = 100 * fully_specific_count / total_designs if total_designs else 0
        print(f"\nDesign-level specificity: {fully_specific_count}/{total_designs} ({designs_pct:.1f}%)")
        print(f"  - {fully_specific_count} designs have all components specific")
        print(f"  - {total_designs - fully_specific_count} designs have at least one cross-reactive component")
        
        # Show where the results were saved
        print(f"\nResults saved to: {results_out_path}")
        print(f"  - File shows all designs with binding site analysis for each component")
        print(f"  - All components are evaluated using appropriate criteria:")
        print(f"    - PCR primers: assessed using detailed binding site criteria")
        print(f"    - QDB probes: assessed using standard identity/coverage criteria")
        print(f"  - Detailed output includes alignment details, Delta G values, and 3' mismatches")
        print(f"  - Sequences are shown for cross-reactive hits to help troubleshoot designs")

if __name__ == "__main__":
    main()
