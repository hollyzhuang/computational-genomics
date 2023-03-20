from collections import Counter
from bwt import *
from compsci260lib import *


'''def main():
    query = "TA"
    reference = "GATTACTA"
    bwt_data = make_all(reference)
    output = find(query, bwt_data)
    print(output)'''


def find(query, bwt_data):
    """Given a query sequence and a series of data structures containing
    various information about the reference text, return a list containing
    all the locations of the query sequence in the reference text.

    Args:
        query (str): query sequence to identify in the reference text
        bwt_data (tuple): BWT data structures containing information about
        the reference text. Refer to the documentation for `make_all` for
        more information.

    Returns:
        (list of ints): of all locations of the query sequence in the
                        reference text
    """

    bwt, suffix_array, ranks, counts = bwt_data

    length = len(bwt)
    results = []
    #
    # YOUR CODE GOES HERE
    #
    # i: the index of the last character of the query
    i = len(query) - 1
    # s: the last character of the query
    s = query[i]
    # sp is the starting point of the last character of the query on the sorted first column
    sp = counts[s] + 1
    # ep is the ending point, which is acquired by adding the total number of appearances of that character to the starting point
    ep = sp + ranks[s][length-1]
    # iterate to the second last character of the query
    i = i - 1
    # while loop ends until the last character is examined or the query is not found at all
    while i >= 0 and sp < ep:
        s = query[i]
        # check the starting point for the next character and locate their position based on ranks
        sp = (counts[s] + 1) + ranks[s][sp - 1]
        ep = (counts[s] + 1) + ranks[s][ep - 1]
        # iterate to the next character
        i = i - 1
    # if the index returned is sp = ep, then it means that the query is not found (e.g.: index 4 to index 4)
    if sp >= ep:
        return "not found"
    else:
        # when found, record the index locations based on suffix array
        pre_results = suffix_array[sp:ep]
        for each in pre_results:
            results.append(each)
    # sort the returned index obtained from the suffix array
    return sorted(results)


# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.

def rank(bwt_seq):
    """Takes as input a string transformed by the BWT. Returns a dictionary
    with characters as keys and lists as values. Each list contains the total
    number of occurrences for the corresponding character up until each
    position in the BWT-transformed string (i.e., its rank).

    For example:
        rank('ACTGA$TA')['$'] --> [0, 0, 0, 0, 0, 1, 1, 1]
        rank('ACTGA$TA')['A'] --> [1, 1, 1, 1, 2, 2, 2, 3]
        rank('ACTGA$TA')['C'] --> [0, 1, 1, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['G'] --> [0, 0, 0, 1, 1, 1, 1, 1]
        rank('ACTGA$TA')['T'] --> [0, 0, 1, 1, 1, 1, 2, 2]
    """
    rank = {}
    characters = set(bwt_seq)
    for character in characters:
        rank[character] = [0]
    rank[bwt_seq[0]] = [1]
    for letter in bwt_seq[1:]:
        for k, v in list(rank.items()):
            v.append(v[-1] + (k == letter))
    return rank


def make_suffix_array(seq):
    """Returns the suffix array of the input string.  NOTE: This is NOT 
    an efficient implementation; it is possible to make a suffix array 
    in much less time asymptotically, but it's a bit more complex.

    For example: print make_suffix_array('GATTACA$') --> [7, 6, 4, 1, 5, 0, 3, 2]
    """
    suffixes = {}
    for x in range(len(seq)):
        suffixes[seq[x:]] = x
    suffix_array = [suffixes[suffix] for suffix in sorted(suffixes.keys())]
    return suffix_array


def count_smaller_chars(seq):
    """Takes as input a string. Returns a dictionary with characters as keys
    and integers as values. The integers track the number of characters in the
    input string which are lexicographically smaller than the corresponding
    character key.

    For example, using an input DNA sequence like 'GATTACA':
        count_smaller_chars('GATTACA')['A'] --> 0
            (A, being lexicographically first in a DNA sequence,
            should always return 0)

        count_smaller_chars('GATTACA')['C'] --> 3
            (C, being second, should return the number of A's, which here is 3)

        count_smaller_chars('GATTACA')['G'] --> 4
            (G, being third, should return the number of A's or C's,
            which here is 4)

        count_smaller_chars('GATTACA')['T'] --> 5
            (T, being fourth, should return the number of A's or C's or G's,
            which here is 5)
    """
    characters = set(seq)
    cntr = Counter(seq)
    total = 0
    counts = {}
    for character in sorted(characters):
        counts[character] = total
        total += cntr[character]
    return counts


def make_all(reference):
    """Takes as input a reference text.

    Returns the data structures necessary to perform efficient exact
    string matching searches.
    """
    counts = count_smaller_chars(reference)
    reference = reference + '$'
    suffix_array = make_suffix_array(reference)
    bwt = forward_bwt(reference)
    ranks = rank(bwt)
    return bwt, suffix_array, ranks, counts


'''if __name__ == '__main__':
    main()'''
