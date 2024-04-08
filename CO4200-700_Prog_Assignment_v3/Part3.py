import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment
from sklearn.cluster import AgglomerativeClustering


# Function to calculate Kimura distance between two alignments
def calculate_kimura_distance(alignment_x, alignment_y):
    positions_scored = 0
    exact_matches = 0

    # Iterate over the alignments
    for a, b in zip(alignment_x, alignment_y):
        # Only consider positions where neither sequence has a gap
        if a != '-' and b != '-':
            positions_scored += 1
            # Count the number of exact matches
            if a == b:
                exact_matches += 1

    # If no positions were scored, return infinity
    if positions_scored == 0:
        return float('inf')  # Avoid division by zero; no valid comparison can be made.

    # Calculate the proportion of exact matches
    S = exact_matches / positions_scored
    # Calculate the proportion of mismatches
    D = 1 - S

    try:
        # Calculate the Kimura distance
        distance = -np.log(1 - D - 0.2 * D**2)
    except ValueError:
        # If D is too close to 1, the log's argument is out of domain
        distance = float('inf')
    
    return distance

# Function to compute alignments for a set of sequences in a file
def compute_alignments(file_path):
    # Open the file and read the lines
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        
    # Get the number of sequences from the first line
    num_sequences = int(lines[0].strip())
    # Get the sequences from the remaining lines
    sequences = [line.strip() for line in lines[1:num_sequences + 1]]

    # Initialize a matrix to store the distances
    distance_matrix = np.zeros((num_sequences, num_sequences))

    
    # Iterate over all pairs of sequences
    for i in range(num_sequences):
        for j in range(num_sequences):
            # Compute the alignment for the pair of sequences
            alignment = Alignment(sequences[i], sequences[j])
            alignment.compute_alignment()
            
            # Get the aligned sequences
            aligned_x, aligned_y = alignment.get_aligned_sequences()
            
            # Calculate the Kimura distance and store it in the matrix
            distance_matrix[i, j] = calculate_kimura_distance(aligned_x, aligned_y)

    # Return the distance matrix
    return distance_matrix


# Function to perform hierarchical clustering on a distance matrix
def cluster(arr):
    # Initialize the clustering model
    model = AgglomerativeClustering(linkage='average', metric='precomputed', distance_threshold = None)
    # Fit the model to the data
    cluster = model.fit(arr)
    # Get the number of objects (sequences)
    n_objects = len(cluster.labels_)
    next_node = n_objects
    output = ""
    
    # Iterate over the merges performed by the clustering
    for i, merge in enumerate(cluster.children_):
        # Add the merge to the output string
        output += f"Align {merge[0]} with {merge[1]} to give {next_node}\n"
        next_node+=1
    # Return the output string
    return output.strip()

# Specify the file path
file_path = "CO4200-700_Prog_Assignment_v3\\sequences\\multiple128.txt"
# Compute the alignments for the sequences in the file
question3i = compute_alignments(file_path)
# Perform hierarchical clustering on the distance matrix
question3ii = cluster(question3i)

# Uncomment the following lines to print out my answers
#print(question3i)
#print(question3ii)