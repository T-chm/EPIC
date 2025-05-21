import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from docker_chat import run_mafft

def read_aligned_sequences(folder_path):
    sequences = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta"):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, "r") as file:
                seq_records = list(SeqIO.parse(file, "fasta"))
                sequences[filename] = seq_records
    return sequences

def sliding_window(seq, window_size, step_size):
    for i in range(0, len(seq) - window_size + 1, step_size):
        window_seq = seq[i:i + window_size]
        if len(window_seq) == window_size and '-' not in window_seq and 'N' not in window_seq and 'n' not in window_seq:
            yield window_seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def design_qdb_assays(sequences, log_file, output_file):
    window_size = 100
    step_size = 100

    with open(log_file, "w") as log, open(output_file, "w") as output:
        for filename, seq_records in sequences.items():
            log.write(f"Processing file: {filename}\n")
            for record in seq_records:
                log.write(f"Analyzing sequence: {record.id}\n")
                for window_seq in sliding_window(str(record.seq), window_size, step_size):
                    target_seq = window_seq.upper()
                    if len(target_seq) != window_size:
                        continue

                    # Split target sequence into two halves
                    half_size = len(target_seq) // 2
                    first_half = target_seq[:half_size]
                    second_half = target_seq[half_size:]

                    # Reverse complement each half
                    reporter_probe = reverse_complement(first_half)
                    capture_probe = reverse_complement(second_half)

                    # Check GC content and melting temperature
                    if reporter_probe and capture_probe:
                        reporter_gc = gc_fraction(reporter_probe) * 100
                        capture_gc = gc_fraction(capture_probe) * 100
                        reporter_tm = mt.Tm_NN(Seq(reporter_probe))
                        capture_tm = mt.Tm_NN(Seq(capture_probe))

                        # Check design criteria
                        if (40 <= reporter_gc <= 60 and 40 <= capture_gc <= 60 and
                            55 <= reporter_tm <= 72 and 55 <= capture_tm <= 72):
                            # Add modifications
                            capture_probe = "/5AmMC6/" + capture_probe
                            reporter_probe = reporter_probe + "/3Cy5Sp/"

                            # Write to output
                            output.write(f"Target: {target_seq}\n")
                            output.write(f"Capture Probe: {capture_probe} (Length: {len(capture_probe)}, GC: {capture_gc:.2f}%, Tm: {capture_tm:.2f}°C)\n")
                            output.write(f"Reporter Probe: {reporter_probe} (Length: {len(reporter_probe)}, GC: {reporter_gc:.2f}%, Tm: {reporter_tm:.2f}°C)\n")
                            output.write(f"Sequence ID: {record.id}\n")
                            output.write(f"Description: {record.description}\n")
                            output.write("\n")

if __name__ == "__main__":
    folder_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "aligned_genome_for_validation_assays")
    log_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qdb_design_log.txt")
    output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qdb_assays_output.txt")

    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    sequences = read_aligned_sequences(folder_path)
    design_qdb_assays(sequences, log_file, output_file)