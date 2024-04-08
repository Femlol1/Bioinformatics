import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment
from Part3 import compute_alignments
from sklearn.cluster import AgglomerativeClustering


# Define a class for nodes in the guide tree
class GuideTreeNode:
    def __init__(self, index, left=None, right=None):
        self.index = index # Index of the sequence this node represents
        self.left = left  # Left child node
        self.right = right # Right child node

# Function to compute progressive alignment
def compute_progressive_alignment(sequences, guide_tree, alignment_instance):
    # If the guide tree node is a leaf node, return the sequence it represents
    if guide_tree.left is None and guide_tree.right is None:
        return [sequences[guide_tree.index]]  # Return a list containing one sequence.

    # Recursively compute alignments for the left and right subtrees
    left_alignment = compute_progressive_alignment(sequences, guide_tree.left, alignment_instance) if guide_tree.left else []
    right_alignment = compute_progressive_alignment(sequences, guide_tree.right, alignment_instance) if guide_tree.right else []

    # Merge the two alignments
    return alignment_instance.computeProfileAlignment(left_alignment, right_alignment)

# Function to build a guide tree from hierarchical clustering
def build_guide_tree(clustering, n_leaves):
    # Initialize the leaf nodes
    nodes = [GuideTreeNode(i) for i in range(n_leaves)]
    next_index = n_leaves

    # Build the internal nodes
    for merge in clustering.children_:
        left = nodes[merge[0]]
        right = nodes[merge[1]]
        nodes.append(GuideTreeNode(next_index, left, right))
        next_index += 1

    return nodes[-1]  # Return the root of the guide tree


# Main function
def main():
    # Specify the file path
    file_path = 'CO4200-700_Prog_Assignment_v3/sequences/multiple128.txt'
    
    # Read the number of sequences and the sequences from the file
    num_sequences = open(file_path, 'r').readlines()[0]
    sequences = [line.strip() for line in open(file_path, 'r').readlines()[1:]]  # Skip the first line

    # Compute the distance matrix
    D = compute_alignments(file_path)
    
    # Perform hierarchical clustering
    model = AgglomerativeClustering(n_clusters=None, linkage='average', metric='precomputed', distance_threshold=0)
    cluster = model.fit(D)

    # Build the guide tree
    guide_tree = build_guide_tree(cluster, len(sequences))

    # Initialize an instance of MultipleAlignment
    alignment_instance = MultipleAlignment()
    
    # Compute the progressive alignment
    multiple_alignment = compute_progressive_alignment(sequences, guide_tree, alignment_instance)
    
    # Print the results
    print(f"Reading input sequences from file {file_path}")
    print(f"Number of sequences {num_sequences}")
    print("Computed alignment:")
    alignment_instance.displayAlignment(multiple_alignment)
    print("Score of alignment: " + str(alignment_instance.scoreMultipleAlignment(multiple_alignment)))


if __name__ == "__main__":
    main()
