import os

from NW_Part1 import Alignment

script_path = os.path.abspath(__file__)

# Get the directory containing the script
script_dir = os.path.dirname(script_path)

# Construct the path to the sequences file
filename = os.path.join(script_dir, 'sequences', 'multiple3.txt')

def nw_align(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """Needleman-Wunsch algorithm for sequence alignment."""
    m, n = len(seq1), len(seq2)
    # Initialize scoring matrix
    score = np.zeros((m+1, n+1), dtype=int)
    # Initialize gap penalties in the first row and column
    for i in range(m+1):
        score[i][0] = i * gap
    for j in range(n+1):
        score[0][j] = j * gap
    # Fill in the scoring matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap
            score[i][j] = max(diag, delete, insert)
    return score[m][n]

def compute_alignment_matrix(filename):
    """Computes the NW alignment score matrix for all sequence pairs."""
    
    # Read sequences from file
    with open(filename, 'r') as f:
        sequences = f.readlines()
    n = len(sequences)
    # Initialize the score matrix
    score_matrix = np.zeros((n, n))
    # Compute scores for all pairs
    for i in range(n):
        for j in range(n):
            score_matrix[i][j] = nw_align(sequences[i], sequences[j])
    return score_matrix

# Example usage
# filename = "./sequences/multiple3.txt"  # Adjust the path as needed
score_matrix = compute_alignment_matrix(filename)
print(f"Computing Alignment for {len(score_matrix)} Sequences in {filename}...")
print(f"Total Comparisons: {len(score_matrix)**2}")
print(score_matrix)
