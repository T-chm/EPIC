import os
import primer3
import multiprocessing
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

def read_qdb_assays(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    assays = []
    assay = {}
    for line in lines:
        line = line.strip()
        if line.startswith("Target: "):
            assay['target'] = line.replace("Target: ", "")
        elif line.startswith("Capture Probe: "):
            assay['capture_probe'] = line.replace("Capture Probe: ", "")
        elif line.startswith("Reporter Probe: "):
            assay['reporter_probe'] = line.replace("Reporter Probe: ", "")
        elif line.startswith("GenBank ID: "):
            assay['genbank_id'] = line.replace("GenBank ID: ", "")
        elif line.startswith("Metadata: "):
            assay['metadata'] = line.replace("Metadata: ", "")
            assays.append(assay)
            assay = {}
    return assays

def design_pcr_primers(assay):
    target_sequence = assay['target']
    seq_args = {
        'SEQUENCE_ID': assay['genbank_id'],
        'SEQUENCE_TEMPLATE': target_sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(target_sequence)]
    }
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 30,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_INTERNAL_MAX_POLY_X': 4,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[80, 100]]
    }
    
    try:
        primers = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as e:
        print(f"Error designing primers for {assay['genbank_id']}: {e}")
        return None
    
    results = []
    for i in range(10):
        try:
            fwd_primer = primers[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_primer = primers[f'PRIMER_RIGHT_{i}_SEQUENCE']
            fwd_tm = primers[f'PRIMER_LEFT_{i}_TM']
            rev_tm = primers[f'PRIMER_RIGHT_{i}_TM']
            fwd_gc = gc_fraction(fwd_primer) * 100
            rev_gc = gc_fraction(rev_primer) * 100
            amplicon_start = primers[f'PRIMER_LEFT_{i}'][0]
            amplicon_end = primers[f'PRIMER_RIGHT_{i}'][0]
            amplicon_seq = target_sequence[amplicon_start:amplicon_end + 1]
            amplicon_rev_comp = str(Seq(amplicon_seq).reverse_complement())
            amplicon_length = len(amplicon_seq)
            results.append({
                'fwd_primer': fwd_primer,
                'rev_primer': rev_primer,
                'fwd_tm': fwd_tm,
                'rev_tm': rev_tm,
                'fwd_gc': fwd_gc,
                'rev_gc': rev_gc,
                'amplicon_seq': amplicon_seq,
                'amplicon_rev_comp': amplicon_rev_comp,
                'amplicon_length': amplicon_length,
                'original_target': target_sequence,
                'metadata': assay['metadata']
            })
        except KeyError:
            break
    return results

def write_primer_to_files(primer, log_file, output_file):
    """Write a single primer result to both log and output files"""
    log_file.write(f"Forward Primer: {primer['fwd_primer']} (Tm: {primer['fwd_tm']}, GC: {primer['fwd_gc']}%)\n")
    log_file.write(f"Reverse Primer: {primer['rev_primer']} (Tm: {primer['rev_tm']}, GC: {primer['rev_gc']}%)\n")
    log_file.write(f"Amplicon Sequence: {primer['amplicon_seq']} (Length: {primer['amplicon_length']})\n")
    log_file.write(f"Amplicon Reverse Complement: {primer['amplicon_rev_comp']}\n")
    log_file.write(f"Original Target: {primer['original_target']}\n")
    log_file.write(f"Metadata: {primer['metadata']}\n\n")
    
    output_file.write(f"Forward Primer: {primer['fwd_primer']} (Tm: {primer['fwd_tm']}, GC: {primer['fwd_gc']}%)\n")
    output_file.write(f"Reverse Primer: {primer['rev_primer']} (Tm: {primer['rev_tm']}, GC: {primer['rev_gc']}%)\n")
    output_file.write(f"Amplicon Sequence: {primer['amplicon_seq']} (Length: {primer['amplicon_length']})\n")
    output_file.write(f"Amplicon Reverse Complement: {primer['amplicon_rev_comp']}\n")
    output_file.write(f"Original Target: {primer['original_target']}\n")
    output_file.write(f"Metadata: {primer['metadata']}\n\n")


def process_and_write_result(assay, log_file_path, output_file_path):
    """Process a single assay and write results directly to files"""
    result = design_pcr_primers(assay)
    
    # Use a file lock to ensure thread-safe file writing
    file_lock = multiprocessing.Lock()
    
    with file_lock:
        with open(log_file_path, 'a') as log_file, open(output_file_path, 'a') as output_file:
            if result is None:
                return
            
            for primer in result:
                write_primer_to_files(primer, log_file, output_file)

def main():
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'qdb_assays_output_all_genomes.txt')
    log_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pcr_assay_log.txt')
    output_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pcr_assays_output.txt')
    
    # Create/clear the output files before starting
    with open(log_file_path, 'w') as f:
        f.write("# PCR Primer Design Log\n\n")
    with open(output_file_path, 'w') as f:
        f.write("# PCR Primer Design Results\n\n")
    
    assays = read_qdb_assays(file_path)
    
    # Determine the number of processes to use
    num_processes = multiprocessing.cpu_count()
    print(f"Using {num_processes} processes for primer design")
    
    # Create a partial function with the file paths
    process_func = partial(process_and_write_result, log_file_path=log_file_path, output_file_path=output_file_path)
    
    # Process each assay and stream results to files using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        # Submit all tasks and get futures
        futures = [executor.submit(process_func, assay) for assay in assays]
        
        # Process results as they complete
        total = len(futures)
        completed = 0
        for future in futures:
            # This will block until the future is done, but we process results as they complete
            future.result()
            completed += 1
            print(f"Progress: {completed}/{total} assays processed ({completed/total*100:.1f}%)")

if __name__ == '__main__':
    main()