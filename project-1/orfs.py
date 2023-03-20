import re

from compsci260lib import get_fasta_dict


def run_orfs():
    """Report the number of ORFs if the minimum ORF length is 10 amino acids,
    40 amino acids, or 60 amino acids.  Also report the average length (in
    amino acids) of the identified ORFs for the three cases.
    """

    #
    # Your code goes here
    #
    # save the file name
    sars_fasta_path = "sars_cov2_wu.fasta"

    # acquire the sars sequence
    sars_dict = get_fasta_dict(sars_fasta_path)
    sars_seq = sars_dict["MN908947.3"].upper()

    # report the number of ORF identified and average ORF length
    min_aalength_10 = 10
    min_aalength_40 = 40
    min_aalength_60 = 60
    orfs_list_10 = find_orfs(sars_seq, min_aalength_10)
    orfs_list_40 = find_orfs(sars_seq, min_aalength_40)
    orfs_list_60 = find_orfs(sars_seq, min_aalength_60)
    summarize_10 = summarize_orfs(orfs_list_10)
    summarize_40 = summarize_orfs(orfs_list_40)
    summarize_60 = summarize_orfs(orfs_list_60)
    print(orfs_list_60)
    print(summarize_60)
    print(f"If the minimum ORF length is {min_aalength_10} amino acids: "
          f"the number of ORF identified is {summarize_10[0]}, with average length in amino acids of {summarize_10[1]}. "
          f"If the minimum ORF length is {min_aalength_40} amino acids: "
          f"the number of ORF identified is {summarize_40[0]}, with average length in amino acids of {summarize_40[1]}. "
          f"If the minimum ORF length is {min_aalength_60} amino acids: "
          f"the number of ORF identified is {summarize_60[0]}, with average length in amino acids of {summarize_60[1]}. ")

    return  # placeholder code, for a valid function


def summarize_orfs(orfs):
    """Summarize ORFs identified from the find_orfs procedure as a count of the
    number of found orfs and the average length of the found ORFs (in amino
    acids)

    Args:
        orfs (list): a list of dictionaries of found ORFs

    Returns:
        tuple: (The number of ORFs found (int), Average ORF length (float))
    """

    #
    # Your code here
    #
    # record the number of ORFs found
    num_orf = len(orfs)
    sum_orf_length = 0
    # calculate the sum of amino acid length for all the ORFs found
    for dict_orf in orfs:
        sum_orf_length += dict_orf["aalength"]
    # get the average by dividing the sum with the number of ORFs found
    ave_orf_length = sum_orf_length/num_orf

    return num_orf, round(ave_orf_length, 2)


def find_orfs(seq, min_length_aa):
    """This is a function for finding sufficiently long ORFs in all reading
    frames in a sequence of DNA or RNA.  By default, the sequence is assumed
    to be single-stranded.

    The function takes as input parameters: a string representing a genomic
    sequence, the minimum length (in amino acids) for an ORF before it will be
    returned (which defaults to 0).

    Args:
        seq (str): a genomic sequence
        min_length_aa (int): minimum length of found ORFs in amino acids

    Returns:
        list: of dictionaries with information on each ORF found.

    Where each ORF found is represented by a dictionary with
    the following keys:
        frame (int): The nucleotide offset in which the ORF was found. (Must be
        0, 1, or 2)
        stop (int): the nucleotide position of the end of the ORF
        start (int): the nucleotide position of the start of the ORF
        stopcodon (str): the nucleotide triplet of the stop codon
        nlength (int): the length (in nucleotides) of the ORF
        strand (str): the strand of the found ORF (Must be "W" or "C")

    A valid return list may look something like this:
    [
        {
            "frame": 0,
            "stop": 13413,
            "aalength": 4382,
            "start": 265,
            "stopcodon": "UAA",
            "nlength": 13149,
            "strand": "W"
        },
        {
            "frame": 0,
            "stop": 27063,
            "aalength": 221,
            "start": 26398,
            "stopcodon": "UAA",
            "nlength": 666,
            "strand": "C"
        }
    ]
    """

    #
    # Your code here
    #
    # create a returned list and store the start and stop codon information
    orf_list = []
    start_codon = "AUG"
    end_c1 = "UAG"
    end_c2 = "UGA"
    end_c3 = "UAA"

    # loop through the whole sequence with increment by 1
    for i in range(len(seq) - 3):
        # find the first AUG
        if seq[i:i + 3] == start_codon:
            # record the start location and the frame
            start = i + 1
            frame = i % 3
            # look for the next (closest) stop codon with increment by 3
            for j in range(1, len(seq) - 3):
                found = False
                if seq[(i + (3 * j)):(i + (3 * j) + 3)] == end_c1 or seq[
                                                                     (i + (3 * j)):(i + (3 * j) + 3)] == end_c2 or seq[(
                        i + (3 * j)):(i + (3 * j) + 3)] == end_c3:
                    # boolean found become true when the next stop codon is found,
                    # record stop, stopcodon, nlength, and aalength values
                    found = True
                    stop = i + (3 * j) + 3
                    stopcodon = seq[i + (3 * j):i + (3 * j) + 3]
                    nlength = 3 * j + 3
                    aalength = j
                    # only include the ORFs that have an amino acid length higher than min length,
                    # create a dictionary, and add to the list
                    if aalength >= min_length_aa:
                        new_dict = {"frame": frame,
                                    "stop": stop,
                                    "aalength": aalength,
                                    "start": start,
                                    "stopcodon": stopcodon,
                                    "nlength": nlength,
                                    "strand": "W"}
                        orf_list.append(new_dict)
                # go back to the first for loop when the closest stop codon is found
                if found:
                    break

    # create a stopcodon dictionary to get the longest ORF if several ORFs share the same
    # stop codon position but have varying lengths
    stopcodon_repeat_dict = {}
    for k in orf_list:
        if k["stop"] not in stopcodon_repeat_dict:
            stopcodon_repeat_dict[k["stop"]] = []
            stopcodon_repeat_dict[k["stop"]].append(k)
        else:
            stopcodon_repeat_dict[k["stop"]].append(k)

    # find the stop codon positions that have more than one ORFs
    # sharing the stop codon with different AUG starting points
    for key in stopcodon_repeat_dict:
        # if there's more than one, only restore the first one — the longest ORF
        if len(stopcodon_repeat_dict[key]) > 1:
            func = lambda x: x["aalength"]
            # reverse order from the longest aalength to the shortest
            stopcodon_repeat_dict[key] = sorted(stopcodon_repeat_dict[key], key=func, reverse=True)
            # store the first element
            stopcodon_repeat_dict[key] = stopcodon_repeat_dict[key][0]
        else:
            # store the first element as well (even tho only have 1 element in the list)
            stopcodon_repeat_dict[key] = stopcodon_repeat_dict[key][0]

    # create a list to store all the ORF dictionaries from the above stop codon dictionary
    new_orf_list = []
    for key in stopcodon_repeat_dict:
        new_orf_list.append(stopcodon_repeat_dict[key])

    return new_orf_list


if __name__ == "__main__":
    """Run run_orfs(). Do not modify this code"""
    run_orfs()
