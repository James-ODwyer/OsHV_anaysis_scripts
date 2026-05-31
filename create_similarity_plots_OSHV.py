# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 16:36:20 2023

@author: odw014
"""

import os
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from tempfile import NamedTemporaryFile
import subprocess


# Specify the path to your desired working directory
new_working_directory = "C:/Users/odw014/OneDrive - CSIRO/Documents/Postdoc bioinformatics/Oyster genome project/final_genomes_maffts_Mauve_trees/All_samples_including_SRAs/whole_genomes"

# Use the os.chdir() function to set the working directory
os.chdir(new_working_directory)




mfa_file = "All_genomes_OSHV_fasta_incl_SRAs.mafft"


alignment = AlignIO.read(mfa_file, "fasta")

# Print basic information about the alignment
print(f"Alignment length: {alignment.get_alignment_length()}")
print(f"Number of sequences: {len(alignment)}")


reference_id = "MF509813.1"
#reference_id = "AY509253.2"
#reference_id = "OM811597.1"
#reference_id = "MW412420.1"
target_ids = ["10-04648-0005", "12-02397-0004", "13-00205-0004", "16-00320", "19-00509", "19-00679-0186", "21-01514", "21-01515", "21-04170", "24-00653-0013", "24-00653-0004", "24-00653-0006","SRR14210303","SRR14210304","SRR14210305"]




def split_alignment(alignment, chunk_size):
    """
    Split a multiple sequence alignment into chunks of a specified size.

    Parameters:
    - alignment: Biopython MultipleSeqAlignment object
    - chunk_size: Size of each chunk

    Returns:
    - List of MultipleSeqAlignment objects representing chunks
    """
    chunks = []
    for start in range(0, len(alignment[0]), chunk_size):
        end = min(start + chunk_size, len(alignment[0]))
        chunk = alignment[:, start:end]
        chunks.append(chunk)
    return chunks

def pairwise_identity(reference, target):
    """
    Calculate pairwise identity between two sequences.

    Parameters:
    - reference: Biopython Seq object
    - target: Biopython Seq object

    Returns:
    - Pairwise identity as a percentage
    """
    alignment_length = min(len(reference), len(target))
    matches = sum(r == t for r, t in zip(reference, target))
    identity = (matches / alignment_length) * 100
    return identity

reference_sequence = [record for record in alignment if record.id == reference_id][0]
target_sequences = [record for record in alignment if record.id in target_ids]

# Split the alignment into 1000 bp chunks
chunk_size = 1000
alignment_chunks = split_alignment(alignment, chunk_size)

# Perform pairwise identity calculations for each chunk
for chunk_num, chunk_alignment in enumerate(alignment_chunks):
    print(f"\nChunk {chunk_num + 1}:")

    # Align the reference and target sequences for the chunk
    reference_chunk_sequence = [record for record in chunk_alignment if record.id == reference_id][0]
    target_chunk_sequences = [record for record in chunk_alignment if record.id in target_ids]


    # Print pairwise identity for each target genome
    for target_id, target_chunk_sequence in zip(target_ids, target_chunk_sequences):
        identity = pairwise_identity(reference_chunk_sequence, target_chunk_sequence)
        print(f"{target_id}: {identity:.2f}%")

"""

Above saves output to terminal

"""



"""

Editing to save to a text file

"""
output_file_path = "Pairwise_identity_OSHV_ref_M509813.txt"
#output_file_path = "Pairwise_identity_OSHV_ref_AY509253.txt"
#output_file_path = "Pairwise_identity_OSHV_ref_OM811597.txt"
#output_file_path = "Pairwise_identity_OSHV_ref_MW412420.txt"




with open(output_file_path, "w") as output_file:
    # Write header
    output_file.write("\t" + "\t".join(map(str, range(1, len(alignment_chunks) + 1))) + "\n")

    # Perform pairwise identity calculations for each chunk
    for target_id in target_ids:
        output_file.write(target_id + "\t")

        # Print pairwise identity for each chunk
        for chunk_num, chunk_alignment in enumerate(alignment_chunks):
            reference_chunk_sequence = [record for record in chunk_alignment if record.id == reference_id][0]
            target_chunk_sequence = [record for record in chunk_alignment if record.id == target_id][0]
            identity = pairwise_identity(reference_chunk_sequence.seq, target_chunk_sequence.seq)
            output_file.write(f"{identity:.2f}\t")

        output_file.write("\n")


