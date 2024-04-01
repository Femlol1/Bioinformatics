import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment
from sklearn.cluster import AgglomerativeClustering


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
        # Log's argument is out of domain, which can happen if D is too close to 1.
        distance = float('inf')
    
    return distance

def compute_alignments(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    num_sequences = int(lines[0].strip())
    sequences = [line.strip() for line in lines[1:num_sequences + 1]]

    #print(f"Computing Kimura Distance for {num_sequences} Sequences in {file_path}...")
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(num_sequences):
            alignment = Alignment(sequences[i], sequences[j])
            alignment.compute_alignment()

            # For demonstration, assume we somehow extract or access these after computation.
            # You will need to replace the next two lines with your method of accessing the sequences.
            aligned_x, aligned_y = alignment.get_aligned_sequences()

            distance_matrix[i, j] = calculate_kimura_distance(aligned_x, aligned_y)

    # print(f"Total Comparisons: {num_sequences ** 2}")
    # print(distance_matrix)
    return distance_matrix


# Modify this path as necessary for your environment.
D = compute_alignments('CO4200-700_Prog_Assignment_v3\sequences\multiple10.txt')

model = AgglomerativeClustering(linkage='average', metric='precomputed', distance_threshold = None)
#cluster the objects based on the distances
cluster = model.fit(D)
n_objects = len(cluster.labels_) # this is just the number of sequences
# cluster.children_ gives the hierarchical clustering
# leaves are the original objects and labelled from 0 to 4 here
# the internal nodes in the cluster are labelled from 5 onwards
next_node = n_objects 
for i, merge in enumerate(cluster.children_):
    print("Align ", merge[0], " with ", merge[1], " to give ", next_node)
    next_node+=1
print("done")
