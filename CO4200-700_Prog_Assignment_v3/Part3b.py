import numpy as np
from multiple_sequence_alignment import MultipleAlignment
from NW_Part1 import Alignment
from Part3 import compute_alignments


def compute_average_distances(D, alignment_order, unaligned_sequences):
    # Calculate the average distance of unaligned sequences to the current alignment
    avg_distances = {}
    for seq_index in unaligned_sequences:
        distances = [D[seq_index, aligned_index] for aligned_index in alignment_order]
        avg_distances[seq_index] = np.mean(distances)
    return avg_distances

def find_initial_pair(D):
    # Find the pair of sequences with the minimum distance
    min_dist = np.inf
    pair = (-1, -1)
    for i in range(D.shape[0]):
        for j in range(i + 1, D.shape[0]):
            if D[i, j] < min_dist:
                min_dist = D[i, j]
                pair = (i, j)
    return pair

def main():
    file_path = 'CO4200-700_Prog_Assignment_v3/sequences/multiple128.txt'
    sequences = [line.strip() for line in open(file_path, 'r').readlines()[1:]]

    D = compute_alignments(file_path)

    # Find the initial closest pair of sequences
    initial_pair = find_initial_pair(D)
    alignment_order = [initial_pair[0], initial_pair[1]]
    unaligned_sequences = set(range(len(sequences))) - set(alignment_order)

    alignment_instance = MultipleAlignment()
    current_alignment = alignment_instance.computeProfileAlignment([sequences[initial_pair[0]]], [sequences[initial_pair[1]]])

    # Iteratively add the remaining sequences based on average distance
    while unaligned_sequences:
        avg_distances = compute_average_distances(D, alignment_order, unaligned_sequences)
        next_sequence = min(avg_distances, key=avg_distances.get)
        current_alignment = alignment_instance.computeProfileAlignment(current_alignment, [sequences[next_sequence]])
        alignment_order.append(next_sequence)
        unaligned_sequences.remove(next_sequence)

    print("Computed alignment:")
    alignment_instance.displayAlignment(current_alignment)
    print("Score of alignment: " + str(alignment_instance.scoreMultipleAlignment(current_alignment)))

if __name__ == "__main__":
    main()
