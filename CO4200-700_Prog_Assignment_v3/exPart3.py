import numpy as np
from NW_Part1 import Alignment


def compute_alignments(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    num_sequences = int(lines[0].strip())  # Number of sequences from the first line
    sequences = [line.strip() for line in lines[1:num_sequences + 1]]  # Sequence data starts from the second line

    if num_sequences != len(sequences):
        print("Error: Number of sequences does not match the expected count.")
        return

    print(f"Computing Alignment for {num_sequences} Sequences in {file_path}...")
    scores_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            alignment = Alignment(sequences[i], sequences[j])
            alignment.compute_alignment()
            scores_matrix[i, j] = alignment.score_alignment()

    print(f"Total Comparisons: {num_sequences ** 2}")
    print(scores_matrix)

# Example usage
compute_alignments('CO4200-700_Prog_Assignment_v3\sequences\multiple128.txt')
