from aligner_helpers import *
from compsci260lib import *
from GlobalAligner import *


def run_global_aligner_plus():
    """Generate the optimal global alignments between:

        2.c)
        atpa_Hs.fasta, atpaEc.fasta

        2.f)
        atpaMm.fasta, atpaHs.fasta
        atpaMm.fasta, atpaEc.fasta

        atpaBs.fasta, atpaHs.fasta
        atpaBs.fasta, atpaEc.fasta

        atpaMm.fasta, atpaBs.fasta

    For each alignment, report the optimal alignment score,
    the top-most alignment, and the bottom-most alignment.
    """

    #
    # YOUR CODE GOES HERE
    #
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
    output = solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty)

    print(output)

    print("The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # 2(c): atpa_Hs.fasta, atpaEc.fasta
    atpa_Hs_dict = get_fasta_dict("atpa_Hs.fasta")
    atpa_Hs_seq = atpa_Hs_dict["sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1"]
    atpa_Ec_dict = get_fasta_dict("atpa_Ec.fasta")
    atpa_Ec_seq = atpa_Ec_dict["sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1"]

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

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Hs_seq, atpa_Ec_seq, subst_dict, gap_penalty)

    print("The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # 2.f)
    # align atpaMm.fasta and atpaHs.fasta
    atpa_Mm_dict = get_fasta_dict("atpa_Mm.fasta")
    atpa_Mm_seq = atpa_Mm_dict["sp|Q03265|ATPA_MOUSE ATP synthase subunit alpha, mitochondrial OS=Mus musculus OX=10090 GN=Atp5f1a PE=1 SV=1"]
    atpa_Hs_dict = get_fasta_dict("atpa_Hs.fasta")
    atpa_Hs_seq = atpa_Hs_dict["sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1"]

    seq_type = validate_sequences(atpa_Mm_seq, atpa_Hs_seq)
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

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Mm_seq, atpa_Hs_seq, subst_dict, gap_penalty)

    print("atpaMm.fasta, atpaHs.fasta - The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # align atpaMm.fasta and atpaEc.fasta
    atpa_Mm_dict = get_fasta_dict("atpa_Mm.fasta")
    atpa_Mm_seq = atpa_Mm_dict["sp|Q03265|ATPA_MOUSE ATP synthase subunit alpha, mitochondrial OS=Mus musculus OX=10090 GN=Atp5f1a PE=1 SV=1"]
    atpa_Ec_dict = get_fasta_dict("atpa_Ec.fasta")
    atpa_Ec_seq = atpa_Ec_dict["sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1"]

    seq_type = validate_sequences(atpa_Mm_seq, atpa_Ec_seq)
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

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Mm_seq, atpa_Ec_seq, subst_dict, gap_penalty)

    print("atpaMm.fasta, atpaEc.fasta - The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # align atpaBs.fasta and atpaHs.fasta
    atpa_Bs_dict = get_fasta_dict("atpa_Bs.fasta")
    atpa_Bs_seq = atpa_Bs_dict["sp|P37808|ATPA_BACSU ATP synthase subunit alpha OS=Bacillus subtilis (strain 168) OX=224308 GN=atpA PE=1 SV=3"]
    atpa_Hs_dict = get_fasta_dict("atpa_Hs.fasta")
    atpa_Hs_seq = atpa_Hs_dict["sp|P0ABB0|ATPA_ECOLI ATP synthase subunit alpha OS=Escherichia coli (strain K12) OX=83333 GN=atpA PE=1 SV=1"]

    seq_type = validate_sequences(atpa_Bs_seq, atpa_Hs_seq)
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

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Bs_seq, atpa_Hs_seq, subst_dict, gap_penalty)

    print("atpaBs.fasta, atpaHs.fasta - The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # align atpaBs.fasta and atpaEc.fasta
    atpa_Bs_dict = get_fasta_dict("atpa_Bs.fasta")
    atpa_Bs_seq = atpa_Bs_dict["sp|P37808|ATPA_BACSU ATP synthase subunit alpha OS=Bacillus subtilis (strain 168) OX=224308 GN=atpA PE=1 SV=3"]
    atpa_Ec_dict = get_fasta_dict("atpa_Ec.fasta")
    atpa_Ec_seq = atpa_Ec_dict["sp|P25705|ATPA_HUMAN ATP synthase subunit alpha, mitochondrial OS=Homo sapiens OX=9606 GN=ATP5F1A PE=1 SV=1"]

    seq_type = validate_sequences(atpa_Bs_seq, atpa_Ec_seq)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit("atpaBs.fasta, atpaEc.fasta - Input sequences are of different types: not both DNA or both protein")

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Bs_seq, atpa_Ec_seq, subst_dict, gap_penalty)

    print("atpaBs.fasta, atpaEc.fasta - The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

    # align atpaMm.fasta and atpaBs.fasta
    atpa_Mm_dict = get_fasta_dict("atpa_Mm.fasta")
    atpa_Mm_seq = atpa_Mm_dict["sp|Q03265|ATPA_MOUSE ATP synthase subunit alpha, mitochondrial OS=Mus musculus OX=10090 GN=Atp5f1a PE=1 SV=1"]
    atpa_Bs_dict = get_fasta_dict("atpa_Bs.fasta")
    atpa_Bs_seq = atpa_Bs_dict["sp|P37808|ATPA_BACSU ATP synthase subunit alpha OS=Bacillus subtilis (strain 168) OX=224308 GN=atpA PE=1 SV=3"]

    seq_type = validate_sequences(atpa_Mm_seq, atpa_Bs_seq)
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

    gap_penalty = 8

    subst_dict = create_subst_matrix_dict(subst_matrix)
    output = solve_global_aligner_plus(atpa_Mm_seq, atpa_Bs_seq, subst_dict, gap_penalty)

    print("atpaMm.fasta, atpaBs.fasta - The optimal alignment score:")
    print(output[0])

    print("The top-most alignment achieving this score:")
    print_alignment(output[1][0], output[1][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)
    print("The bottom-most alignment achieving this score:")
    print_alignment(output[2][0], output[2][1], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

def solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    filling in the DP table, and returning the final value and alignments.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

    Returns a tuple of:
        (the optimal alignment score as an int,
         the top-most alignment achieving this score as a tuple of strings,
         the bottom-most alignment achieving this score as a tuple of strings)

        Example output:

            (6, ("AT-AGG", "ATCCGG"), ("ATA-GG", "ATCCGG"))

    Note: If you do the extra challenge to report all optimal alignments,
    you can lengthen the size of the return tuple, but ensure that its second
    element (i.e., the first alignment) remains the top-most alignment, while
    the last element is the bottom-most alignment. e.g.

        (optimal score, (top-most alignment sequences), ...,
         (bottom-most alignment sequences))
    """

    #
    # YOUR CODE GOES HERE
    #
    # Initialize the DP table's data structure
    # as a list of lists of ints

    dp_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    updated_dp_table = global_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # Compute the score of the optimal global alignment
    max_value = updated_dp_table[len(dp_table) - 1][len(dp_table[0]) - 1]

    # initialize the alignment sequences with their length in blanks
    top_seq_1 = ""
    top_seq_2 = ""
    bottom_seq_1 = ""
    bottom_seq_2 = ""

    # start the traceback from the lower right corner of the dp table
    index_i = len(dp_table) - 1
    index_j = len(dp_table[0]) - 1

    # top-most alignment: up - diagonal - left
    while index_i > 0 and index_j > 0:
        string_match = seq1[index_i - 1] + seq2[index_j - 1]
        # store the penalized scores for each of the type
        last_type_1 = updated_dp_table[index_i - 1][index_j - 1] + subst_dict[string_match]
        last_type_2 = updated_dp_table[index_i - 1][index_j] - gap_penalty
        last_type_3 = updated_dp_table[index_i][index_j - 1] - gap_penalty

        # check whether it equals to type 2 and then type 1 and then type 3 because we want the top-most sequence
        # if one condition is met, then proceed to the next cell with the while loop
        if last_type_2 == updated_dp_table[index_i][index_j]:
            # update sequence alignment and the indexes of cell on the dp_table based on the last alignment type
            top_seq_1 = seq1[index_i - 1] + top_seq_1
            top_seq_2 = "-" + top_seq_2
            index_i = index_i - 1
            index_j = index_j
        elif last_type_1 == updated_dp_table[index_i][index_j]:
            top_seq_1 = seq1[index_i - 1] + top_seq_1
            top_seq_2 = seq2[index_j - 1] + top_seq_2
            index_i = index_i - 1
            index_j = index_j - 1
        elif last_type_3 == updated_dp_table[index_i][index_j]:
            top_seq_1 = "-" + top_seq_1
            top_seq_2 = seq2[index_j - 1] + top_seq_2
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

    # bottom-most alignment: left - diagonal - up
    # start the traceback from the lower right corner of the dp table
    index_i = len(dp_table) - 1
    index_j = len(dp_table[0]) - 1

    # check whether it equals to type 3 and then type 1 and then type 2 because we want the buttom-most sequence
    # if one condition is met, then proceed to the next cell with the while loop
    while index_i > 0 and index_j > 0:
        string_match = seq1[index_i - 1] + seq2[index_j - 1]
        last_type_1 = updated_dp_table[index_i - 1][index_j - 1] + subst_dict[string_match]
        last_type_2 = updated_dp_table[index_i - 1][index_j] - gap_penalty
        last_type_3 = updated_dp_table[index_i][index_j - 1] - gap_penalty

        # update sequence alignment and the indexes of cell on the dp_table based on the last alignment type
        if last_type_3 == updated_dp_table[index_i][index_j]:
            bottom_seq_1 = "-" + bottom_seq_1
            bottom_seq_2 = seq2[index_j - 1] + bottom_seq_2
            index_i = index_i
            index_j = index_j - 1
        elif last_type_1 == updated_dp_table[index_i][index_j]:
            bottom_seq_1 = seq1[index_i - 1] + bottom_seq_1
            bottom_seq_2 = seq2[index_j - 1] + bottom_seq_2
            index_i = index_i - 1
            index_j = index_j - 1
        elif last_type_2 == updated_dp_table[index_i][index_j]:
            bottom_seq_1 = seq1[index_i - 1] + bottom_seq_1
            bottom_seq_2 = "-" + bottom_seq_2
            index_i = index_i - 1
            index_j = index_j

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

    return max_value, (top_seq_1, top_seq_2), (bottom_seq_1, bottom_seq_2)


# new helper function to just return the dp_table for traceback instead of optimal
# score found at the lower right corner of the dp table
def global_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table):
    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)  # so I is 1 more than m
    J = len(dp_table[0])  # so J is 1 more than n

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
    for i in range(1, I):
        for j in range(1, J):
            string_match = seq1[i - 1] + seq2[j - 1]
            type_1 = dp_table[i - 1][j - 1] + subst_dict[string_match]
            type_2 = dp_table[i - 1][j] - gap_penalty
            type_3 = dp_table[i][j - 1] - gap_penalty
            max_value = max(type_1, type_2, type_3)
            dp_table[i][j] = max_value

    # return the dp_table for traceback:
    return dp_table


if __name__ == "__main__":
    run_global_aligner_plus()
