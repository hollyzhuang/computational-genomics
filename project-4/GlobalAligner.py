from aligner_helpers import *
from compsci260lib import *


def run_global_aligner():
    """
    Given two sequences of either DNA or amino acids, initialize the
    appropriate substitution matrix, run the global aligner and report the
    optimal alignment score.
    """

    seq1 = "GAATCGGA"
    seq2 = "TAGTA"
    match = 2
    mismatch = -1
    gap_penalty = 2

    seq_type = validate_sequences(seq1, seq2)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit("Input sequences are of different types: not both DNA or both protein")

    # Obtain a dictionary of scores for aligning a pair of characters
    subst_dict = create_subst_matrix_dict(subst_matrix)

    # store the score
    optimal_score = solve_global_aligner(seq1, seq2, subst_dict, gap_penalty)

    print(optimal_score)

    # Report the optimal global alignment score
    #
    # YOUR CODE GOES HERE
    #


def solve_global_aligner(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and return the optimal alignment score

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): penalty for a column containing a gap char (g: use a
            positive value because this value will be subtracted)

    Returns:
        (int) the optimal alignment score
    """

    # Initialize the DP table's data structure
    # as a list of lists of ints
    dp_table = [[0] * (len(seq2)+1) for _ in range(len(seq1)+1)]

    # Compute the score of the optimal global alignment
    max_value = global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # Display the dp table
    display_dp_table(seq1, seq2, dp_table)

    return max_value


def global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table):
    """A dynamic programming algorithm that takes two sequences and returns the
    score of the optimal alignment.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """

    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)      # so I is 1 more than m
    J = len(dp_table[0])   # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    #
    # YOUR CODE GOES HERE
    #
    for i in range(1, I):
        dp_table[i][0] = dp_table[i][0] - (i * gap_penalty)
    for j in range(1, J):
        dp_table[0][j] = dp_table[0][j] - (j * gap_penalty)

    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    #
    # YOUR CODE GOES HERE
    #
    # fill out the table rows by rows based on the max value from left, diag, and up
    for i in range(1, I):
        for j in range(1, J):
            string_match = seq1[i-1] + seq2[j-1]
            type_1 = dp_table[i-1][j-1] + subst_dict[string_match]
            type_2 = dp_table[i-1][j] - gap_penalty
            type_3 = dp_table[i][j-1] - gap_penalty
            max_value = max(type_1, type_2, type_3)
            # store the max value
            dp_table[i][j] = max_value

    # The optimal score is found at the lower right corner of the dp table:
    return dp_table[I-1][J-1]


def compute_global_aligner_score_linspace(seq1, seq2, subst_dict, gap_penalty):
    """Extra Challenge: compute the score (not the actual alignment) of the
       best global alignment in O(mn) time but only O(min(m,n)) space

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """

    seqshort, seqlong = (seq1, seq2) if len(seq1) < len(seq2) else (seq2, seq1)

    #
    # YOUR CODE GOES HERE
    #
    return -1


if __name__ == "__main__":
    run_global_aligner()
