import numpy as np
from multiple_sequence_alignment import MultipleAlignment
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


# Here is a sample output of working with the NW algorithm from Part1
x = "GSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKL"
y = "NNPELQAHAGKVFKLVYEAAIQLQVTGVVVTDATLKNLGSVHVSKG"

# Here we are calling the alignment method
a = Alignment(x, y) # instantiates the alignment object, providing sequences X and Y
a.compute_alignment() # run to do the actual computation
print("Computed alignment:")
a.display_alignment() # this can be displayed to help visualise what is going on

# a.score_alignment provides a float score for X and Y alignment
print("Score of alignment: " + str(a.score_alignment()))

compute_alignments('CO4200-700_Prog_Assignment_v3\sequences\multiple3.txt')
