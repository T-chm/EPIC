import os
from Bio import AlignIO
from io import StringIO
from docker_chat import run_mafft  # Import the pre-defined MAFFT function

# Define the folder containing the FASTA files
input_folder = 'viral_genomes'
output_folder = 'aligned_viral_genomes'

# Ensure the output directory exists


os.makedirs(output_folder, exist_ok=True)

# Process each FASTA file


for filename in os.listdir(input_folder):


   if filename.endswith('.fasta'):
       filepath = os.path.join(input_folder, filename)
       
        # Read sequences from the FASTA file
       sequences = {}


       with open(filepath, 'r') as file:
           for record in SeqIO.parse(file, 'fasta'):
               sequences[record.id] = str(record.seq)

        # Align sequences using MAFFT


       alignment_result = run_mafft(sequences, output_format="clustal")

        # Parse the alignment result


       alignment_io = StringIO(alignment_result)


       alignment = AlignIO.read(alignment_io, "clustal")

        # Write the alignment to a new file
       output_filepath = os.path.join(output_folder, f"aligned_{filename}")
       with open(output_filepath, 'w') as output_file:
           AlignIO.write(alignment, output_file, "clustal")

print(f"Aligned sequences from {filename} and saved to {output_filepath}")
