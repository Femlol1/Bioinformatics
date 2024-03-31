import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


def read_sequences(filename):
    with open(filename, 'r') as file:
        sequences = file.read().split()
    return sequences

def nw_alignment_score(seq1, seq2, matrix=matlist.blosum62, gap_open=-10, gap_extend=-0.5):
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    # Return the score of the top alignment
    return alignments[0].score

def compute_alignment_matrix(sequences):
    n = len(sequences)
    score_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                # Aligning the sequence with itself, take the full length match score
                score = nw_alignment_score(sequences[i], sequences[j])
            else:
                score = nw_alignment_score(sequences[i], sequences[j])
                score_matrix[j][i] = score
            score_matrix[i][j] = score
    return score_matrix

# Example usage
filename = "./sequences/multiple3.txt"
print(f"Computing Alignment for sequences in {filename}...")
sequences = read_sequences(filename)
score_matrix = compute_alignment_matrix(sequences)
n = len(sequences)
print(f"Total Comparisons: {n**2}")
print(score_matrix)
