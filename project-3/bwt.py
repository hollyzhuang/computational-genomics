from compsci260lib import *

# Note that the '$' character will be used to designate the end of a given
# string.


def solve_bwt():
    """Load the mystery.txt file and decode its contents using the
    reverse BWT transformation. Report the decoded contents."""

    # Example sequences for forward_bwt() and reverse_bwt();
    # you may change them to different sequences to test your code.
    seq1 = 'GATTACA$'
    output1 = forward_bwt(seq1)
    seq2 = 'ACTGA$TA'
    output2 = reverse_bwt(seq2)
    print(output1)
    print(output2)

    # report question (d):
    question_d_sequence = 'CGGACTAACGGACTAACGGACTAACGGACTAC$'
    question_d_output = forward_bwt(question_d_sequence)
    print("The permuted sequence output is: " + question_d_output)

    """
    seq1 = 'GATTACA$'
    seq2 = 'ACTGA$TA'
    forward_bwt(seq1)
    reverse_bwt(seq2)
    """

    # Code to open the file mystery.txt in the correct encoding
    # across platforms, and read its contents into a variable.

    with open('mystery.txt', 'r', encoding='UTF-8') as f:
         mystery_seq = f.read()

    #
    # YOUR CODE GOES HERE
    #
    mystery_output = reverse_bwt(mystery_seq)
    # replace underscore with space
    mystery_output = mystery_output.replace('_', ' ')
    print(mystery_output)


def forward_bwt(seq):
    """forward_bwt(seq) takes as input a string containing the EOF character to
    which the BWT must be applied. The method should then return the result of
    the BWT on the input string.

    For example:
        forward_bwt('GATTACA$') --> 'ACTGA$TA'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """

    #
    # YOUR CODE GOES HERE
    #
    # create an empty suffix array
    suffix_array = []
    # put every suffix of the original sequence in the list
    for i in range(len(seq)):
        suffix_array.append(seq[i:])
    sorted_array = suffix_array
    # sort the array based on the first letter
    sorted_array.sort()
    # in the sorted array, fill up the rest of the sequence after the dollar sign so that they are all rotational
    for i in range(len(sorted_array)):
        length = len(sorted_array[i])
        # based on the different length of the sequence in the sorted array,
        # the rest that needs to be filled up is dependent on the length of the sorted one
        sorted_array[i] += seq[:len(seq)-length]
    transformed = ""
    # get the last column
    for sequence in sorted_array:
        transformed += sequence[-1]
    return transformed

def reverse_bwt(seq):
    """reverse_bwt(seq) takes as input a string containing the EOF character to
    which the reverse of the BWT must be applied. The method should then return
    the result of the reversal on the input string.

    For example:
        reverse_bwt('ACTGA$TA') --> 'GATTACA$'

    Args:
        seq (str): input string with an EOF character

    Returns:
        (str): the transformed string
    """

    #
    # YOUR CODE GOES HERE
    #
    # create a sequence list and a sort list
    seq_list = []
    sort_list = []
    for n in range(len(seq)):
        seq_list.append(seq[n])
    for n in range(len(seq)):
        sort_list.append(seq[n])
    # sort the sequence, which is the first column of the bwt structure
    sort_list.sort()

    # appending
    append_list = []
    for k in range(len(seq) - 1):
        # every time a new append happens, the append-list will be cleared to empty
        append_list.clear()
        # loop through every row
        for i in range(len(seq)):
            # combined is the result of appending the input sequence list in front of the sorted list for each character
            combined = seq_list[i] + sort_list[i]
            # for every character's combination, add it to the append_list
            append_list.append(combined)
        # clear the sorted list
        sort_list.clear()
        # copy the whole append_list to sort_list and then sort the sort_list
        for element in append_list:
            sort_list.append(element)
        sort_list.sort()

    # looping through the list to find the original sequence
    for element in sort_list:
        # the sequence with the dollar sign at the end is the original sequence
        if element[-1] == "$":
            return element


if __name__ == '__main__':
    """Run solve_bwt(). Do not modify this code"""
    solve_bwt()
