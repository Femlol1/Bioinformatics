import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment
from Part3 import compute_alignments
from sklearn.cluster import AgglomerativeClustering


class GuideTreeNode:
    def __init__(self, index, left=None, right=None):
        self.index = index
        self.left = left
        self.right = right

def compute_progressive_alignment(sequences, guide_tree, alignment_instance):
    if guide_tree.left is None and guide_tree.right is None:
        return [sequences[guide_tree.index]]  # Return a list containing one sequence.

    left_alignment = compute_progressive_alignment(sequences, guide_tree.left, alignment_instance) if guide_tree.left else []
    right_alignment = compute_progressive_alignment(sequences, guide_tree.right, alignment_instance) if guide_tree.right else []

    # Merge the two alignments
    return alignment_instance.computeProfileAlignment(left_alignment, right_alignment)

def build_guide_tree(clustering, n_leaves):
    nodes = [GuideTreeNode(i) for i in range(n_leaves)]
    next_index = n_leaves

    for merge in clustering.children_:
        left = nodes[merge[0]]
        right = nodes[merge[1]]
        nodes.append(GuideTreeNode(next_index, left, right))
        next_index += 1

    return nodes[-1]  # Return the root of the guide tree

def main():
    file_path = 'CO4200-700_Prog_Assignment_v3/sequences/multiple10.txt'
    sequences = [line.strip() for line in open(file_path, 'r').readlines()[1:]]  # Skip the first line

    D = compute_alignments(file_path)
    model = AgglomerativeClustering(n_clusters=None, linkage='average', metric='precomputed', distance_threshold=0)
    cluster = model.fit(D)

    guide_tree = build_guide_tree(cluster, len(sequences))

    alignment_instance = MultipleAlignment()
    multiple_alignment = compute_progressive_alignment(sequences, guide_tree, alignment_instance)
    
    print("Computed alignment:")
    alignment_instance.displayAlignment(multiple_alignment)
    print("Score of alignment: " + str(alignment_instance.scoreMultipleAlignment(multiple_alignment)))

if __name__ == "__main__":
    main()
