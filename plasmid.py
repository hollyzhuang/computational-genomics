import re

from compsci260lib import get_fasta_dict, reverse_complement, translate


def run_plasmid():
    """Function to
    1) Run the simple_assembler function to assemble the plasmid from a .fasta of reads
    2) Run the find_orfs_circular_double_stranded function to find ORFs in the assembled
       circular plasmid (considering both Watson and Crick strands)
    3) Report the number of ORFs found that are of >= 50 amino acids long, coding fractions
    4) Translate the longest ORF to its amino acid sequence and print it out.
    """
    # Load the plasmid reads from fasta, and ensure they are in proper order to be assembled
    #
    # Your code here
    #
    plasmid_dict = get_fasta_dict("plasmid.fasta")
    reads = []
    # append the first element in the list to be the read whose key is "start"
    for key in plasmid_dict.keys():
        if key == "start":
            reads.append(plasmid_dict[key])
    # append the rest of the reads to the list who are not "start"
    for key in plasmid_dict.keys():
        if key != "start":
            reads.append(plasmid_dict[key])

    # Assemble the reads into a single-stranded linearized version of the plasmid:
    plasmid_dna = simple_assembler(reads)
    print(plasmid_dna)
    # Report the sequence of the single-stranded linearized version of the plasmid:
    print(f"The complete assembled sequence is 5' -- {plasmid_dna[0:30]} ... {plasmid_dna[8520:]} --3'")

    # Report the length of the assembled plasmid
    #
    # Your code here
    #
    assembled_dna_length = len(plasmid_dna)
    print(assembled_dna_length)

    # Search for ORFs in the reconstructed plasmid DNA and report your results
    #
    # Your code here
    #
    plasmid_orfs = find_orfs_circular_double_stranded(plasmid_dna.upper(), 50)
    print(plasmid_orfs)
    print(len(plasmid_orfs))

    # report the output of the found ORFs
    print(
        f"ORFs that have a minimum length of 50 amino acids include one sequence on the {plasmid_orfs[0]['strand']} strand that starts from {plasmid_orfs[0]['start']} and ends at {plasmid_orfs[0]['stop']} with amino acid length of {plasmid_orfs[0]['aalength']}, "
        f"one sequence on the {plasmid_orfs[1]['strand']} strand that starts from {plasmid_orfs[1]['start']} and ends at {plasmid_orfs[1]['stop']} with amino acid length of {plasmid_orfs[1]['aalength']}, "
        f"one sequence on the {plasmid_orfs[2]['strand']} strand that starts from {plasmid_orfs[2]['start']} and ends at {plasmid_orfs[2]['stop']} with amino acid length of {plasmid_orfs[2]['aalength']}, "
        f"one sequence on the {plasmid_orfs[3]['strand']} strand that starts from {plasmid_orfs[3]['start']} and ends at {plasmid_orfs[3]['stop']} with amino acid length of {plasmid_orfs[3]['aalength']}, "
        f"one sequence on the {plasmid_orfs[4]['strand']} strand that starts from {plasmid_orfs[4]['start']} and ends at {plasmid_orfs[4]['stop']} with amino acid length of {plasmid_orfs[4]['aalength']}, "
        f"one sequence on the {plasmid_orfs[5]['strand']} strand that starts from {plasmid_orfs[5]['start']} and ends at {plasmid_orfs[5]['stop']} with amino acid length of {plasmid_orfs[5]['aalength']}, "
        f"one sequence on the {plasmid_orfs[6]['strand']} strand that starts from {plasmid_orfs[6]['start']} and ends at {plasmid_orfs[6]['stop']} with amino acid length of {plasmid_orfs[6]['aalength']}, "
        f"one sequence on the {plasmid_orfs[7]['strand']} strand that starts from {plasmid_orfs[7]['start']} and ends at {plasmid_orfs[7]['stop']} with amino acid length of {plasmid_orfs[7]['aalength']}, "
        f"one sequence on the {plasmid_orfs[8]['strand']} strand that starts from {plasmid_orfs[8]['start']} and ends at {plasmid_orfs[8]['stop']} with amino acid length of {plasmid_orfs[8]['aalength']}, "
        f"one sequence on the {plasmid_orfs[9]['strand']} strand that starts from {plasmid_orfs[9]['start']} and ends at {plasmid_orfs[9]['stop']} with amino acid length of {plasmid_orfs[9]['aalength']}, "
        f"one sequence on the {plasmid_orfs[10]['strand']} strand that starts from {plasmid_orfs[10]['start']} and ends at {plasmid_orfs[10]['stop']} with amino acid length of {plasmid_orfs[10]['aalength']}, "
        f"one sequence on the {plasmid_orfs[11]['strand']} strand that starts from {plasmid_orfs[11]['start']} and ends at {plasmid_orfs[11]['stop']} with amino acid length of {plasmid_orfs[11]['aalength']}, "
        f"one sequence on the {plasmid_orfs[12]['strand']} strand that starts from {plasmid_orfs[12]['start']} and ends at {plasmid_orfs[12]['stop']} with amino acid length of {plasmid_orfs[12]['aalength']}, "
        f"one sequence on the {plasmid_orfs[13]['strand']} strand that starts from {plasmid_orfs[13]['start']} and ends at {plasmid_orfs[13]['stop']} with amino acid length of {plasmid_orfs[13]['aalength']}, "
        f"one sequence on the {plasmid_orfs[14]['strand']} strand that starts from {plasmid_orfs[14]['start']} and ends at {plasmid_orfs[14]['stop']} with amino acid length of {plasmid_orfs[14]['aalength']}, "
        f"one sequence on the {plasmid_orfs[15]['strand']} strand that starts from {plasmid_orfs[15]['start']} and ends at {plasmid_orfs[15]['stop']} with amino acid length of {plasmid_orfs[15]['aalength']}, "
        f"one sequence on the {plasmid_orfs[16]['strand']} strand that starts from {plasmid_orfs[16]['start']} and ends at {plasmid_orfs[16]['stop']} with amino acid length of {plasmid_orfs[16]['aalength']}, "
        f"one sequence on the {plasmid_orfs[17]['strand']} strand that starts from {plasmid_orfs[17]['start']} and ends at {plasmid_orfs[17]['stop']} with amino acid length of {plasmid_orfs[17]['aalength']}.")

    # Find the longest ORF and translate
    #
    # Your code here
    #
    longest = 0
    for orf in plasmid_orfs: # loop through the list of dictionary to find the ODF with the longest length
        length = orf["aalength"]
        if length > longest:
            longest = length # store the current highest ORF
            return_orf = orf
    print(return_orf)

    # print the longest ORF sequence of the plasmid:
    first_half = plasmid_dna.upper()[8045:] # because the ORF contains sequence at where the plasmid is "cut", print the first half of the sequence
    second_half = plasmid_dna.upper()[:2021-3] # print the second half starting from index 0
    longest_orf = first_half + second_half
    print(longest_orf)

    # translate the longest ORF sequence to amino acid sequence:
    aa_sequence = translate(longest_orf)
    print(aa_sequence)

    # Compute the fraction of the genome that is coding
    fraction_list = []
    # make a string with all zeros
    for i in range(8421):
        fraction_list.append(0)
    # loop through watson strand
    for i in range(0, 8):
        start_i = plasmid_orfs[i]['start'] - 1
        stop_i = plasmid_orfs[i]['stop'] - 3 - 1
        for j in range(start_i, stop_i+1):
            if fraction_list[j] is 0:
                fraction_list[j] = 1
            else:
                fraction_list[j] = 1
    # loop through crick strand
    for i in range(9, 18):
        start_i = plasmid_orfs[i]['start'] - 1
        stop_i = plasmid_orfs[i]['stop'] - 3 - 1
        for j in range(len(fraction_list) - stop_i+1, len(fraction_list) - start_i):
            if fraction_list[j] is 0:
                fraction_list[j] = 1
            else:
                fraction_list[j] = 1
    count = 0
    for i in range(len(fraction_list)):
        if fraction_list[i] == 1:
            count += 1

    fraction = count/8421
    print(fraction)

def simple_assembler(reads):
    """Given a list of reads, as described in the problem, return an assembled
    DNA sequence, as a string. For consistency, use the first entry in the
    fragment reads list as the starting position of the returned sequence.

    For example, if we were to take in a list of three reads, 31 nucleotides
    long each. The last 15 nucleotides of each read would overlap with one
    other read, and the assembled sequence would be 48 nucleotides long with
    the sequence starting with the beginning of the first read.

    Args:
        reads (list): list of sequence strings as reads

    Returns:
         str: an assembled genomic sequence as a string starting with the first
              read in `reads`
    """

    #
    # Your code here
    #
    # put the first sequence in the returned sequence string
    return_sequence = reads[0]
    # algorithm that continuously decreases the size of the input list until it has only one element left in the list
    while len(reads) > 1:
        length_read = len(reads[0])
        # record the overlap sequence
        overlap = reads[0][length_read - 15:]
        # delete the first read in the list
        reads = reads[1:]
        # loop through the rest of the reads in the list and find the one with overlapping region
        for i in range(len(reads)):
            if overlap in reads[i]:
                next_sequence = reads[i]
                # add the read with the overlapping region to the returned string
                return_sequence += next_sequence[15:]
                # swap the read with overlapping region with the first read in the list
                reads[0], reads[i] = reads[i], reads[0]
    # delete the last 15 nucleotides of the last read bc it overlaps with the first 15 nucleotides of the first read
    return_sequence = return_sequence[:len(return_sequence)-15]

    return return_sequence

def find_orfs_circular_double_stranded(seq, min_length_aa=0):
    """This is a function for finding sufficiently long ORFs in all reading
    frames in a sequence of double-stranded circular DNA.

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
        frame (int): the reading frame (or nucleotide offset) in which
            the ORF was found (must be 0, 1, or 2)
        start (int): the nucleotide position of the start of the ORF
            (relative to 5' end of the found ORF's strand)
        stop (int): the nucleotide position of the end of the ORF
            (relative to 5' end of the found ORF's strand)
        stopcodon (str): the nucleotide triplet of the stop codon
        nlength (int): the length (in nucleotides) of the ORF
        aalength (int): the length (in amino acids) of the translated protein
        strand (str): the strand of the found ORF (must be "W" or "C")

    A valid return list may look something like this:
    [
        {
            "frame": 0,
            "stop": 13413,
            "aalength": 4382,
            "start": 265,
            "stopcodon": "TAA",
            "nlength": 13149,
            "strand": "W"
        },
        {
            "frame": 0,
            "stop": 27063,
            "aalength": 221,
            "start": 26398,
            "stopcodon": "TAG",
            "nlength": 666,
            "strand": "C"
        }
    ]
    """

    #
    # Your code here
    #
    # list for watson
    orf_list = []
    # list for crick
    orf_list_complement = []
    # double the size of watson to account for circular plasmid
    seq_double = seq + seq
    # get the reverse complement of watson strand
    seq_complement = reverse_complement(seq)
    # double the size of crick to account for circular plasmid
    seq_complement_double = seq_complement + seq_complement

    start_codon = "ATG"
    end_c1 = "TAG"
    end_c2 = "TGA"
    end_c3 = "TAA"

    # loop through the whole sequence with increment by 1
    for i in range(len(seq) - 3):
        # find the first AUG
        if seq[i:i + 3] == start_codon:
            # record the start location and the frame
            start = i + 1
            frame = i % 3
            # look for the next (closest) stop codon with increment by 3
            for j in range(1, len(seq_double) - 3 - i):
                found = False
                if seq_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c1 or seq_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c2 or seq_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c3:
                    # boolean found become true when the next stop codon is found,
                    # record stop, stopcodon, nlength, and aalength values
                    found = True
                    stop = i + (3 * j) + 3
                    # if the stop index is bigger than the length of the original sequence before
                    # it gets doubled, return the original index
                    if stop > len(seq):
                        stop = stop - len(seq)
                    stopcodon = seq_double[i + (3 * j):i + (3 * j) + 3]
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

    # the same code above but searching ORFs in the Crick(complement) strand
    for i in range(len(seq_complement) - 3):
        if seq_complement[i:i + 3] == start_codon:
            start = i + 1
            frame = i % 3
            for j in range(1, len(seq_complement_double) - 3 - i):
                found = False
                if seq_complement_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c1 or seq_complement_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c2 or seq_complement_double[(i + (3 * j)):(i + (3 * j) + 3)] == end_c3:
                    found = True
                    stop = i + (3 * j) + 3
                    if stop > len(seq_complement):
                        stop = stop - len(seq_complement)
                    stopcodon = seq_complement_double[i + (3 * j):i + (3 * j) + 3]
                    nlength = 3 * j + 3
                    aalength = j
                    if aalength >= min_length_aa:
                        new_dict = {"frame": frame,
                                    "stop": stop,
                                    "aalength": aalength,
                                    "start": start,
                                    "stopcodon": stopcodon,
                                    "nlength": nlength,
                                    "strand": "C"} # change W to C because this is for searching in the C strand
                        orf_list_complement.append(new_dict)
                if found:
                    break

    # same code above but for Crick
    stopcodon_repeat_dict_complement = {}
    for k in orf_list_complement:
        if k["stop"] not in stopcodon_repeat_dict_complement:
            stopcodon_repeat_dict_complement[k["stop"]] = []
            stopcodon_repeat_dict_complement[k["stop"]].append(k)
        else:
            stopcodon_repeat_dict_complement[k["stop"]].append(k)

    # same code for Crick
    for key in stopcodon_repeat_dict_complement:
        if len(stopcodon_repeat_dict_complement[key]) > 1:
            func = lambda x: x["aalength"]
            stopcodon_repeat_dict_complement[key] = sorted(stopcodon_repeat_dict_complement[key], key=func, reverse=True)
            stopcodon_repeat_dict_complement[key] = stopcodon_repeat_dict_complement[key][0]
        else:
            stopcodon_repeat_dict_complement[key] = stopcodon_repeat_dict_complement[key][0]

    # create a list to store all the dictionaries from Watson and Crick
    new_orf_list = []
    for key in stopcodon_repeat_dict:
        new_orf_list.append(stopcodon_repeat_dict[key])
    for key in stopcodon_repeat_dict_complement:
        new_orf_list.append(stopcodon_repeat_dict_complement[key])

    return new_orf_list


if __name__ == "__main__":
    run_plasmid()
