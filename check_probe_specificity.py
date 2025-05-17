#!/usr/bin/env python3
"""
Script to check the specificity of probe sequences against a local BLAST database.
This script parses probe sequences from rep_qdb_assays.txt and runs BLAST searches
against a local virus database to evaluate their specificity.
"""

import os
import re
import subprocess
import csv
import tempfile
from collections import defaultdict

# Configuration
INPUT_FILE = "rep_qdb_assays.txt"
BLAST_DB_PATH = "local_virus_blast_db/ref_viruses_rep_genomes"
BLAST_OUTPUT_DIR = "blast_results"
E_VALUE_THRESHOLD = 1.0  # Maximum E-value to consider as a significant hit
IDENTITY_THRESHOLD = 80  # Minimum percent identity to consider as a significant hit
QUERY_COVERAGE_THRESHOLD = 70  # Minimum query coverage to consider as a significant hit


def parse_probe_file(file_path):
    """Parse the probe file and extract probe sequences."""
    probes = []
    current_probe = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith("Target:"):
            if current_probe:
                probes.append(current_probe)
                current_probe = {}
            
            current_probe["target_seq"] = line.replace("Target:", "").strip()
        
        elif line.startswith("Capture Probe:"):
            # Extract sequence between /5AmMC6/ and space
            match = re.search(r'/5AmMC6/([ATCG]+)', line)
            if match:
                current_probe["capture_probe"] = match.group(1)
            
            # Extract GC content and Tm if available
            gc_match = re.search(r'GC: ([\d.]+)%', line)
            tm_match = re.search(r'Tm: ([\d.]+)', line)
            
            if gc_match:
                current_probe["capture_gc"] = float(gc_match.group(1))
            if tm_match:
                current_probe["capture_tm"] = float(tm_match.group(1))
        
        elif line.startswith("Reporter Probe:"):
            # Extract sequence between beginning and /3Cy5Sp/
            match = re.search(r'^Reporter Probe: ([ATCG]+)/3Cy5Sp/', line)
            if match:
                current_probe["reporter_probe"] = match.group(1)
            
            # Extract GC content and Tm if available
            gc_match = re.search(r'GC: ([\d.]+)%', line)
            tm_match = re.search(r'Tm: ([\d.]+)', line)
            
            if gc_match:
                current_probe["reporter_gc"] = float(gc_match.group(1))
            if tm_match:
                current_probe["reporter_tm"] = float(tm_match.group(1))
        
        elif line.startswith("GenBank ID:"):
            current_probe["genbank_id"] = line.replace("GenBank ID:", "").strip()
        
        elif line.startswith("Metadata:"):
            current_probe["metadata"] = line.replace("Metadata:", "").strip()
        
        i += 1
    
    # Add the last probe if exists
    if current_probe:
        probes.append(current_probe)
    
    return probes


def run_blast_search(sequence, seq_type, probe_info):
    """Run a BLAST search for a given sequence against the local database."""
    # Create a temporary file for the sequence
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as temp:
        temp_path = temp.name
        genbank_id = probe_info.get("genbank_id", "Unknown")
        temp.write(f">{genbank_id}_{seq_type}\n{sequence}\n")
    
    # Ensure output directory exists
    os.makedirs(BLAST_OUTPUT_DIR, exist_ok=True)
    
    # Output file for this blast search
    output_file = os.path.join(BLAST_OUTPUT_DIR, f"{genbank_id}_{seq_type}_blast.tsv")
    
    # Run BLAST command
    blast_cmd = [
        "blastn",
        "-query", temp_path,
        "-db", BLAST_DB_PATH,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle",
        "-evalue", str(E_VALUE_THRESHOLD),
        "-out", output_file,
        "-max_target_seqs", "10"  # Limit to top 10 hits for brevity
    ]
    
    try:
        subprocess.run(blast_cmd, check=True)
        print(f"BLAST search completed for {genbank_id} {seq_type} probe")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        return None
    finally:
        # Clean up the temporary file
        os.unlink(temp_path)
    
    # Parse BLAST results
    results = []
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        with open(output_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 13:  # Ensure we have at least the required columns
                    hit = {
                        'qseqid': parts[0],
                        'sseqid': parts[1],
                        'pident': float(parts[2]),
                        'length': int(parts[3]),
                        'qstart': int(parts[6]),
                        'qend': int(parts[7]),
                        'evalue': float(parts[10]),
                        'bitscore': float(parts[11]),
                        'qcovs': float(parts[12]),
                        'stitle': parts[13] if len(parts) > 13 else ''
                    }
                    results.append(hit)
    
    return results


def analyze_specificity(probe_results):
    """
    Analyze BLAST results to determine probe specificity.
    A probe is considered specific if it only has significant hits against its intended target.
    """
    specificity = {}
    
    for probe_id, results in probe_results.items():
        genbank_id, probe_type = probe_id.split('_', 1)
        
        # Skip if no results
        if not results:
            specificity[probe_id] = {
                'is_specific': True,
                'reason': 'No significant hits found',
                'hits_count': 0
            }
            continue
        
        # Count significant hits that meet our thresholds
        significant_hits = [
            hit for hit in results
            if (hit['pident'] >= IDENTITY_THRESHOLD and
                hit['qcovs'] >= QUERY_COVERAGE_THRESHOLD and
                hit['evalue'] <= E_VALUE_THRESHOLD)
        ]
        
        # Check if the probe hits its intended target
        target_hits = [
            hit for hit in significant_hits
            if genbank_id.split(".")[0] in hit['sseqid']
        ]
        
        # Analyze specificity
        if len(significant_hits) == 0:
            specificity[probe_id] = {
                'is_specific': True,
                'reason': 'No significant hits found',
                'hits_count': 0
            }
        elif len(target_hits) > 0 and len(significant_hits) == len(target_hits):
            specificity[probe_id] = {
                'is_specific': True,
                'reason': 'Only hits intended target',
                'hits_count': len(significant_hits),
                'top_hit': significant_hits[0] if significant_hits else None
            }
        else:
            off_target_hits = [hit for hit in significant_hits if hit not in target_hits]
            specificity[probe_id] = {
                'is_specific': False,
                'reason': f'Hits {len(off_target_hits)} off-target sequences',
                'hits_count': len(significant_hits),
                'off_target_count': len(off_target_hits),
                'top_hit': significant_hits[0] if significant_hits else None,
                'top_off_target': off_target_hits[0] if off_target_hits else None
            }
    
    return specificity


def write_summary_report(probes, probe_results, specificity_analysis):
    """Write a summary report of probe specificity."""
    output_file = os.path.join(BLAST_OUTPUT_DIR, "probe_specificity_summary.csv")
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow([
            'GenBank ID', 'Probe Type', 'Specific', 'Reason',
            'Total Significant Hits', 'Top Hit Identity', 'Top Hit Coverage',
            'Top Hit E-value', 'Top Hit Accession', 'Top Hit Description'
        ])
        
        # Write data for each probe
        for probe in probes:
            genbank_id = probe['genbank_id']
            
            for probe_type in ['capture_probe', 'reporter_probe']:
                if probe_type in probe:
                    probe_id = f"{genbank_id}_{probe_type}"
                    analysis = specificity_analysis.get(probe_id, {
                        'is_specific': False,
                        'reason': 'No analysis available',
                        'hits_count': 0
                    })
                    
                    top_hit = analysis.get('top_hit', {})
                    
                    writer.writerow([
                        genbank_id,
                        probe_type,
                        'Yes' if analysis.get('is_specific', False) else 'No',
                        analysis.get('reason', ''),
                        analysis.get('hits_count', 0),
                        top_hit.get('pident', '') if top_hit else '',
                        top_hit.get('qcovs', '') if top_hit else '',
                        top_hit.get('evalue', '') if top_hit else '',
                        top_hit.get('sseqid', '') if top_hit else '',
                        top_hit.get('stitle', '') if top_hit else ''
                    ])
    
    print(f"Summary report written to {output_file}")


def main():
    # Create output directory if it doesn't exist
    os.makedirs(BLAST_OUTPUT_DIR, exist_ok=True)
    
    # Parse the probe file
    print(f"Parsing probe file {INPUT_FILE}...")
    probes = parse_probe_file(INPUT_FILE)
    print(f"Found {len(probes)} probe sets")
    
    # Run BLAST searches for each probe
    print("Running BLAST searches...")
    probe_results = defaultdict(list)
    
    for i, probe in enumerate(probes):
        print(f"Processing probe set {i+1}/{len(probes)} ({probe.get('genbank_id', 'Unknown')})")
        
        # Process capture probe
        if 'capture_probe' in probe:
            capture_results = run_blast_search(
                probe['capture_probe'],
                'capture_probe',
                probe
            )
            if capture_results:
                probe_id = f"{probe['genbank_id']}_capture_probe"
                probe_results[probe_id] = capture_results
        
        # Process reporter probe
        if 'reporter_probe' in probe:
            reporter_results = run_blast_search(
                probe['reporter_probe'],
                'reporter_probe',
                probe
            )
            if reporter_results:
                probe_id = f"{probe['genbank_id']}_reporter_probe"
                probe_results[probe_id] = reporter_results
    
    # Analyze specificity
    print("Analyzing probe specificity...")
    specificity_analysis = analyze_specificity(probe_results)
    
    # Write summary report
    print("Writing summary report...")
    write_summary_report(probes, probe_results, specificity_analysis)
    
    print("Analysis complete!")


if __name__ == "__main__":
    main()
