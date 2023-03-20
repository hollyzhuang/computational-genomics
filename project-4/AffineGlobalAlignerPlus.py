from compsci260lib import *
from aligner_helpers import *


def run_ag_aligner_plus():
    """Align atpa_Hs.fasta and atpaEc.fasta and report the optimal
    alignment score, the top-most alignment, and the bottom-most
    alignment.
    """

    #
    # YOUR CODE GOES HERE
    #
    seq1 = "GAATCGGA"
    seq2 = "TAGTA"
    match = 2
    mismatch = -1
    gap_penalty = 1
    affine_penalty = 11

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

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty)
    print(output)

    # atpa Hs.fasta and atpa Ec.fasta
    atpa_Hs_dict = get_fasta_dict("atpa_Hs.fasta")
    atpa_Hs_seq = atpa_Hs_dict[
        "sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1"]
    atpa_Ec_dict = get_fasta_dict("atpa_Ec.fasta")
    atpa_Ec_seq = atpa_Ec_dict[
        "sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1"]

    seq_type = validate_sequences(atpa_Hs_seq, atpa_Ec_seq)
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

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_ag_aligner_plus(atpa_Hs_seq, atpa_Ec_seq, subst_dict, gap_penalty, affine_penalty)

    print("The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)


def solve_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty):
    """The procedure for collecting the inputs, running the aligner,
    and returning the score and optimal alignments.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary
            representation of the substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character);
            this value should be positive because we will subtract it
        affine_penalty (int): affine penalty; as a positive integer

    Returns a tuple of:
        (the optimal alignment score as an int,
         the top-most alignment achieving this score as a tuple of
         strings, the bottom-most alignment achieving this score as a
         tuple of strings)

        Example output:
            (6, ("AT-AGG", "ATCCGG"), ("ATA-GG", "ATCCGG"))
    """

    #
    # YOUR CODE GOES HERE
    #
    # initialize the three dp_tables with zeros
    dp_table_1 = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    dp_table_2 = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    dp_table_3 = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # put the dp_tables into the format of V(k, i, j)
    # so we know the first index of dp_table will indicate table D, E, F (0, 1, 2)
    dp_table = []
    dp_table.append(dp_table_1)
    dp_table.append(dp_table_2)
    dp_table.append(dp_table_3)

    # updated_table returns the dp_table
    updated_table = ag_aligner_output_dp_table(seq1, seq2, subst_dict, gap_penalty, affine_penalty, dp_table)

    I = len(dp_table[0])
    J = len(dp_table[0][0])

    # calculate the max score by going to the bottom right of each table and see which table has the highest score
    max_score = float('-inf')
    last_type = -1
    for k in range(3):
        score = updated_table[k][I - 1][J - 1]
        if score > max_score:
            max_score = score
            last_type = k

    # initialize the sequences with their length in blanks
    top_seq_1 = ""
    top_seq_2 = ""
    bottom_seq_1 = ""
    bottom_seq_2 = ""

    # top-most alignment: up - diagonal - left
    # start the traceback from the lower right corner of the dp table with the last alignment with the type scoring the highest
    index_i = I - 1
    index_j = J - 1
    while index_i > 0 and index_j > 0:
        # check which type is the last type first
        if last_type == 0:
            # then record the corresponding sequences
            top_seq_1 = seq1[index_i - 1] + top_seq_1
            top_seq_2 = seq2[index_j - 1] + top_seq_2
            string_match = seq1[index_i - 1] + seq2[index_j - 1]
            last_type_1 = updated_table[0][index_i - 1][index_j - 1] + subst_dict[string_match]
            last_type_2 = updated_table[1][index_i - 1][index_j - 1] + subst_dict[string_match]
            last_type_3 = updated_table[2][index_i - 1][index_j - 1] + subst_dict[string_match]
            if last_type_2 == updated_table[0][index_i][index_j]:
                last_type = 1
            elif last_type_1 == updated_table[0][index_i][index_j]:
                last_type = 0
            elif last_type_3 == updated_table[0][index_i][index_j]:
                last_type = 2
            index_i = index_i - 1
            index_j = index_j - 1
        elif last_type == 1:
            top_seq_1 = seq1[index_i - 1] + top_seq_1
            top_seq_2 = "-" + top_seq_2
            last_type_1 = updated_table[0][index_i - 1][index_j] - affine_penalty - gap_penalty
            last_type_2 = updated_table[1][index_i - 1][index_j] - gap_penalty
            last_type_3 = updated_table[2][index_i - 1][index_j] - affine_penalty - gap_penalty
            if last_type_2 == updated_table[1][index_i][index_j]:
                last_type = 1
            elif last_type_1 == updated_table[1][index_i][index_j]:
                last_type = 0
            elif last_type_3 == updated_table[1][index_i][index_j]:
                last_type = 2
            index_i = index_i - 1
            index_j = index_j
        elif last_type == 2:
            top_seq_1 = "-" + top_seq_1
            top_seq_2 = seq2[index_j - 1] + top_seq_2
            last_type_1 = updated_table[0][index_i][index_j - 1] - affine_penalty - gap_penalty
            last_type_2 = updated_table[1][index_i][index_j - 1] - affine_penalty - gap_penalty
            last_type_3 = updated_table[2][index_i][index_j - 1] - gap_penalty
            if last_type_2 == updated_table[2][index_i][index_j]:
                last_type = 1
            elif last_type_1 == updated_table[2][index_i][index_j]:
                last_type = 0
            elif last_type_3 == updated_table[2][index_i][index_j]:
                last_type = 2
            index_i = index_i
            index_j = index_j - 1
    # when any of index_i or index_j hits zero, we hit the first row or column of the dp_table
    # this means that we know for the longer sequence, the rest of the alignment type will be type 2/type 3
    while index_i > 0:
        top_seq_1 = seq1[index_i - 1] + top_seq_1
        top_seq_2 = "-" + top_seq_2
        index_i -= 1
    while index_j > 0:
        top_seq_1 = "-" + top_seq_1
        top_seq_2 = seq2[index_j - 1] + top_seq_2
        index_j -= 1

    # restart the last_type with bottom-most because last_type could be altered
    max_score = float('-inf')
    last_type = -1
    for k in range(3):
        score = updated_table[k][I - 1][J - 1]
        if score > max_score:
            max_score = score
            last_type = k

    # bottom-most alignment: up - diagonal - left
    index_i = I - 1
    index_j = J - 1
    while index_i > 0 and index_j > 0:
        # check which type is the last type first
        if last_type == 0:
            # then record the corresponsing sequences
            bottom_seq_1 = seq1[index_i - 1] + bottom_seq_1
            bottom_seq_2 = seq2[index_j - 1] + bottom_seq_2
            string_match = seq1[index_i - 1] + seq2[index_j - 1]
            last_type_1 = updated_table[0][index_i - 1][index_j - 1] + subst_dict[string_match]
            last_type_2 = updated_table[1][index_i - 1][index_j - 1] + subst_dict[string_match]
            last_type_3 = updated_table[2][index_i - 1][index_j - 1] + subst_dict[string_match]
            if last_type_3 == updated_table[0][index_i][index_j]:
                last_type = 2
            elif last_type_1 == updated_table[0][index_i][index_j]:
                last_type = 0
            elif last_type_2 == updated_table[0][index_i][index_j]:
                last_type = 1
            index_i = index_i - 1
            index_j = index_j - 1
        elif last_type == 1:
            bottom_seq_1 = seq1[index_i - 1] + bottom_seq_1
            bottom_seq_2 = "-" + bottom_seq_2
            last_type_1 = updated_table[0][index_i - 1][index_j] - affine_penalty - gap_penalty
            last_type_2 = updated_table[1][index_i - 1][index_j] - gap_penalty
            last_type_3 = updated_table[2][index_i - 1][index_j] - affine_penalty - gap_penalty
            if last_type_3 == updated_table[1][index_i][index_j]:
                last_type = 2
            elif last_type_1 == updated_table[1][index_i][index_j]:
                last_type = 0
            elif last_type_2 == updated_table[1][index_i][index_j]:
                last_type = 1
            index_i = index_i - 1
            index_j = index_j
        elif last_type == 2:
            bottom_seq_1 = "-" + bottom_seq_1
            bottom_seq_2 = seq2[index_j - 1] + bottom_seq_2
            last_type_1 = updated_table[0][index_i][index_j - 1] - affine_penalty - gap_penalty
            last_type_2 = updated_table[1][index_i][index_j - 1] - affine_penalty - gap_penalty
            last_type_3 = updated_table[2][index_i][index_j - 1] - gap_penalty
            if last_type_3 == updated_table[2][index_i][index_j]:
                last_type = 2
            elif last_type_1 == updated_table[2][index_i][index_j]:
                last_type = 0
            elif last_type_2 == updated_table[2][index_i][index_j]:
                last_type = 1
            index_i = index_i
            index_j = index_j - 1
    # when any of index_i or index_j hits zero, we hit the first row or column of the dp_table
    # this means that we know for the longer sequence, the rest of the alignment type will be type 2/type 3
    while index_i > 0:
        bottom_seq_1 = seq1[index_i - 1] + bottom_seq_1
        bottom_seq_2 = "-" + bottom_seq_2
        index_i -= 1
    while index_j > 0:
        bottom_seq_1 = "-" + bottom_seq_1
        bottom_seq_2 = seq2[index_j - 1] + bottom_seq_2
        index_j -= 1

    return max_score, (top_seq_1, top_seq_2), (bottom_seq_1, bottom_seq_2)

# function to return the dp_table:
def ag_aligner_output_dp_table(seq1, seq2, subst_dict, gap_penalty, affine_penalty, dp_table):
    I = len(dp_table[0])  # so I is 1 more than m
    J = len(dp_table[0][0])  # so J is 1 more than n

    # Initialize the dp tables with solutions to base cases using affine gap penalty
    # V(k, i, j) = the score of the optimal alignment of the i-prefix of seq1 to the j-prefix of seq2 such that
    # the last column is type k (k only has three possibilities so the three-dimensional matrix can be considered
    # as 3 two-dimensional matrices

    # dp_table_1 (k = 0)
    dp_table[0][0][0] = 0
    for i in range(1, I):
        dp_table[0][i][0] = float('-inf')
    for j in range(1, J):
        dp_table[0][0][j] = float('-inf')

    # dp_table_2 (k = 1)
    for j in range(J):
        dp_table[1][0][j] = float('-inf')
    for i in range(1, I):
        dp_table[1][i][0] = dp_table[1][i][0] - (i * gap_penalty) - affine_penalty

    # dp_table_3 (k = 2)
    for i in range(I):
        dp_table[2][i][0] = float('-inf')
    for j in range(1, J):
        dp_table[2][0][j] = dp_table[2][0][j] - (j * gap_penalty) - affine_penalty

    # Compute the scores for the rest of the matrix
    for i in range(1, I):
        for j in range(1, J):
            # dp_table_1 - type 1
            t1_type_1 = dp_table[0][i - 1][j - 1]
            t1_type_2 = dp_table[1][i - 1][j - 1]
            t1_type_3 = dp_table[2][i - 1][j - 1]
            t1_max_value = max(t1_type_1, t1_type_2, t1_type_3)
            t1_string_match = seq1[i - 1] + seq2[j - 1]
            dp_table[0][i][j] = subst_dict[t1_string_match] + t1_max_value

            # dp_table_2 - type 2
            t2_type_1 = dp_table[0][i - 1][j] - affine_penalty
            t2_type_2 = dp_table[1][i - 1][j]
            t2_type_3 = dp_table[2][i - 1][j] - affine_penalty
            t2_max_value = max(t2_type_1, t2_type_2, t2_type_3)
            dp_table[1][i][j] = t2_max_value - gap_penalty

            # dp_table_3 - type 3
            t3_type_1 = dp_table[0][i][j - 1] - affine_penalty
            t3_type_2 = dp_table[1][i][j - 1] - affine_penalty
            t3_type_3 = dp_table[2][i][j - 1]
            t3_max_value = max(t3_type_1, t3_type_2, t3_type_3)
            dp_table[2][i][j] = t3_max_value - gap_penalty

    return dp_table

if __name__ == "__main__":
    run_ag_aligner_plus()
