import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment


# Define the function to compute alignments
def compute_alignments(file_path):
    # Open the file in read mode
    with open(file_path, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

    # Gets the number of sequences from the first line of the file
    num_sequences = int(lines[0].strip())
    # Get the sequence data starting from the second line
    sequences = [line.strip() for line in lines[1:num_sequences + 1]]  # Sequence data starts from the second line

    # Check if the number of sequences matches the expected count
    if num_sequences != len(sequences):
        print("Error: Number of sequences does not match the expected count.")
        return
    
    # Print the number of sequences and the file path
    # print(f"Computing Alignment for {num_sequences} Sequences in {file_path}...")
    # Initialize a matrix to store the scores
    scores_matrix = np.zeros((num_sequences, num_sequences))

    # Loop through each sequence
    for i in range(num_sequences):
        for j in range(num_sequences):
            # Compute the alignment for each pair of sequences
            alignment = Alignment(sequences[i], sequences[j])
            alignment.compute_alignment()
            # Store the score in the matrix
            scores_matrix[i, j] = alignment.score_alignment()

    # Returns the total number of comparisons and the scores matrix
    
    
    total_comparisons = num_sequences ** 2

    return total_comparisons, scores_matrix,num_sequences,file_path
    

# Call the function with the path to the file
total_comparisons, scores_matrix,num_sequences,file_path = compute_alignments('CO4200-700_Prog_Assignment_v3\\sequences\\multiple32.txt')
# I just wanted to get the exact same result as shown in sample output part2 is at least close enough

relative_path = file_path.split("CO4200-700_Prog_Assignment_v3\\", 1)[1]

print(f"Computing Alignment for {num_sequences} Sequences in {relative_path}...")
print(f"Total Comparisons: {total_comparisons}")
print(scores_matrix)