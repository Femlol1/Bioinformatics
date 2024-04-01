import numpy as np
from NW_Part1 import Alignment


def calculate_kimura_distance(alignment_x, alignment_y):
    positions_scored = 0
    exact_matches = 0

    for a, b in zip(alignment_x, alignment_y):
        if a != '-' and b != '-':
            positions_scored += 1
            if a == b:
                exact_matches += 1

    if positions_scored == 0:
        return float('inf')  # Avoid division by zero; no valid comparison can be made.

    S = exact_matches / positions_scored
    D = 1 - S

    try:
        distance = -np.log(1 - D - 0.2 * D**2)
    except ValueError:
        distance = float('inf')  # Assign infinity if the log's argument is out of domain.
    
    return distance

def compute_alignments(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    num_sequences = int(lines[0].strip())
    sequences = [line.strip() for line in lines[1:num_sequences + 1]]

    print(f"Computing Kimura Distance for {num_sequences} Sequences in {file_path}...")
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            alignment = Alignment(sequences[i], sequences[j])
            alignment.compute_alignment()
            aligned_x, aligned_y = alignment.get_aligned_sequences()  # Assuming this method exists.
            distance_matrix[i, j] = calculate_kimura_distance(aligned_x, aligned_y)

    print(f"Total Comparisons: {num_sequences ** 2}")
    print(distance_matrix)

# Example usage, assuming the file path is compatible with your OS and environment.
compute_alignments('CO4200-700_Prog_Assignment_v3\sequences\multiple3.txt')
