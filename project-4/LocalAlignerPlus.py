from compsci260lib import *
from aligner_helpers import *


def run_local_aligner_plus():
    """Locally align O18381.fasta and P63015.fasta and report the
    optimal score and optimal local alignment information.
    """
    #
    # YOUR CODE GOES HERE
    #
    # print report for the two protein sequences P63015 and O18381
    O18381_dict = get_fasta_dict("O18381.fasta")
    O18381_seq = O18381_dict[
        "sp|O18381|PAX6_DROME Paired box protein Pax-6 OS=Drosophila melanogaster OX=7227 GN=ey PE=1 SV=3"]
    P63015_dict = get_fasta_dict("P63015.fasta")
    P63015_seq = P63015_dict[
        "sp|P63015|PAX6_MOUSE Paired box protein Pax-6 OS=Mus musculus OX=10090 GN=Pax6 PE=1 SV=1"]

    match = 2
    mismatch = -1
    gap_penalty = 8
    seq_type = validate_sequences(O18381_seq, P63015_seq)
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
    output = solve_local_aligner_plus(O18381_seq, P63015_seq, subst_dict, gap_penalty)
    print(output)

    print("The optimal alignment score: " + str(output[0]))

    print("The total number of locations in the table achieving this optimal score is " + str(len(output[1])))

    print("The locations in the table of each of these optima are:")
    print(output[1])

    print("Indices of the characters of the first aligned sequence (inclusive) are: " + str(output[2][0]))
    print("Indices of the characters of the second aligned sequence (inclusive) are: " + str(output[2][1]))

    print("One of the optimal alignments:")
    print_alignment(output[2][2], output[2][3], seq1_start=1, seq2_start=1, name1="Seq 1", name2="Seq 2",
                    nt_per_line=60)

def solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The procedure for collecting the inputs, running the aligner,
    and returning the score and optimal alignment(s).

    Note that for each returned local alignment, starting positions
    also need to be returned. These are the positions of the first
    character in each aligned sequence relative to the original
    sequence.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary
            representation of the substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character);
            this value should be positive because we will subtract it

    A max score may be in multiple locations, so return the optimal
    score, the locations of all the maxima, and any one optimal
    alignment as a tuple.

    Returns a tuple of:
        (the optimal alignment score as an int,

         the locations of the maxima in the dp table as a list of
         tuples. these positions will include the offset of the
         initialized penalty row and column, so that location (i,j)
         refers to the i-prefix of X and the j-prefix of Y, just as in
         lecture,

         tuple for an optimal alignment)

        The alignment will be in the form:

              (tuple of indices of the characters of the first aligned
               sequence used in the alignment),

              (tuple of indices of the characters of the second aligned
               sequence used in the alignment),

              the first aligned sequence as a string,

              the second aligned sequence as a string)

        As an example with the sequences:

            Sequence 1: TAG
            Sequence 2: TAAGATAAG

        A possible return may be:

            (11, # the optimal score

             # the two maximal locations in the dp table
             [(3, 4), (3, 9)],

             # one possible alignment:
             ((1, 3), # the nt positions mapping TA-G from TAG
              (1, 4), # the nt positions mapping TAAG from TAAGATAAG

              "TA-G", "TAAG") # the sequences of the alignment
            )

        Corresponding to two maxima with an optimal
        alignment score of 11.
    """
    #
    # YOUR CODE GOES HERE
    #
    dp_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    updated_dp_table = local_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # updated_dp_table is a tuple storing the dp_table and the max score
    dp_table = updated_dp_table[0]
    opt_align_score = updated_dp_table[1]

    I = len(dp_table)
    J = len(dp_table[0])

    # store the locations of the max alignment scores for local alignment
    maxima_locations = []
    for i in range(I):
        for j in range(J):
            if dp_table[i][j] == opt_align_score:
                location = (i, j)
                maxima_locations.append(location)

    # need to do some traceback here
    # initialize the sequences with their length in blanks
    top_seq_1 = ""
    top_seq_2 = ""
    
    # start the traceback from the indexes of the first-found location with the highest local alignment score
    # before we are only ask to return any alignment, which one doesn't matter
    index_i = maxima_locations[0][0]
    index_j = maxima_locations[0][1]

    # top-most alignment: up - diagonal - left
    # will always do trace-back for the top-most sequence since we're asked to return any one optimal alignment
    # check whether it equals to type 2 and then type 1 and then type 3 because we want the top-most sequence
    # if one condition is met, then proceed to the next cell with the while loop
    while index_i > 0 and index_j > 0:
        string_match = seq1[index_i - 1] + seq2[index_j - 1]
        last_type_1 = dp_table[index_i - 1][index_j - 1] + subst_dict[string_match]
        last_type_2 = dp_table[index_i - 1][index_j] - gap_penalty
        last_type_3 = dp_table[index_i][index_j - 1] - gap_penalty

        # if all left, diag, and up cells of the index cell are zero, we want to proceed the traceback because we want
        # the longer sequence alignment
        # because the way to go is arbitrary - find any sequence alignment
        # we choose to always go up and then diag and then left, so we will always get the top-most alignment
        # type 2 - type 1 - type 3
        if dp_table[index_i][index_j] == 0 and dp_table[index_i - 1][index_j - 1] == 0 and dp_table[index_i - 1][index_j] == 0 and dp_table[index_i][index_j - 1] == 0:
            if last_type_2 == 0:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = "-" + top_seq_2
                index_i = index_i - 1
                index_j = index_j
            elif last_type_1 == 0:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i - 1
                index_j = index_j - 1
            elif last_type_3 == 0:
                top_seq_1 = "-" + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i
                index_j = index_j - 1
        elif dp_table[index_i][index_j] == 0:
            # if the index cell is zero, and all other three cells (left, diag, up) after penalization are negative values,
            # we stop at this index cell of value 0
            if last_type_2 < 0 and last_type_1 < 0 and last_type_3 < 0:
                return_index_i = index_i
                return_index_j = index_j
                break;
            # if the index cell is zero, and there are other cells after penalization that's not zero, then we don't break
            # the while loop and proceed to the next because this zero must come from some other values
            # again, in sequence type 2 -1 - 3 because we want the top-most sequence alignment
            if last_type_2 > 0:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = "-" + top_seq_2
                index_i = index_i - 1
                index_j = index_j
            elif last_type_1 > 0:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i - 1
                index_j = index_j - 1
            elif last_type_3 > 0:
                top_seq_1 = "-" + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i
                index_j = index_j - 1
        # similar process as the GlobalALignerPlus for traceback when the index cell is not 0.
        # in the order of type 2 -1 - 3 because of looking for the top-most alignment
        elif dp_table[index_i][index_j] != 0:
            if last_type_2 == dp_table[index_i][index_j]:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = "-" + top_seq_2
                index_i = index_i - 1
                index_j = index_j
            elif last_type_1 == dp_table[index_i][index_j]:
                top_seq_1 = seq1[index_i - 1] + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i - 1
                index_j = index_j - 1
            elif last_type_3 == dp_table[index_i][index_j]:
                top_seq_1 = "-" + top_seq_1
                top_seq_2 = seq2[index_j - 1] + top_seq_2
                index_i = index_i
                index_j = index_j - 1

        # whenever the traceback hits the first column or the first row, break the while loop
        if index_i == 0 or index_j == 0:
            return_index_i = index_i
            return_index_j = index_j
            break;

    return opt_align_score, maxima_locations, ((return_index_i + 1, maxima_locations[0][0]), (return_index_j + 1, maxima_locations[0][1]), top_seq_1, top_seq_2)

    # example format of return tuple:
    # return (opt_align_score, maxima_locations, ((-1, -1), (-1, -1), "", "")

# helper function to return
# 1. the whole dp_table
# 2. the max score found from the whole dp table
def local_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table):
    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)  # so I is 1 more than m
    J = len(dp_table[0])  # so J is 1 more than n

    #
    # YOUR CODE GOES HERE
    #
    # fill in the table
    return_max = 0
    for i in range(1, I):
        for j in range(1, J):
            string_match = seq1[i - 1] + seq2[j - 1]
            type_1 = dp_table[i - 1][j - 1] + subst_dict[string_match]
            type_2 = dp_table[i - 1][j] - gap_penalty
            type_3 = dp_table[i][j - 1] - gap_penalty
            # added a condition for type 4 = 0 for local alignment because we want to start over whenever
            # all three other cells are negative after penalization
            type_4 = 0
            max_value = max(type_1, type_2, type_3, type_4)
            dp_table[i][j] = max_value
            # to store the max value found from the local alignment dp_table
            if max_value > return_max:
                return_max = max_value

    # return the dp_table for traceback:
    return dp_table, return_max

if __name__ == "__main__":
    run_local_aligner_plus()
