from compsci260lib import *
from pprint import pprint

def run_assembly_tester():
    """Read in the read and contig files, verify the reads in each contig, and estimate
    the gap length between the contigs.
    """
    reads_file = "paired.reads.fasta"
    contigs_file = "contigs.fasta"

    # Load the two fasta files
    reads = get_fasta_dict(reads_file)
    contigs = get_fasta_dict(contigs_file)

    # Determine how the reads map to the contigs
    contig_reads = find_contig_reads(reads, contigs)
    #print(contig_reads)

    # Report:
    #   - Any reads that are not consistent with the overall assembly based
    #     on the checks done by find_contig_reads()
    #   - The read pairs where both reads match to (distinct) contigs in supercontig1,
    #   and their estimated gap lengths
    #
    # YOUR CODE HERE
    #
    # - Any reads that are not consistent with the overall assembly based
    #     on the checks done by find_contig_reads()
    seq_not_belong = []
    for seq in contig_reads:
        if contig_reads[seq]['contig_a'] is None and contig_reads[seq]['contig_b'] is None:
            failed_reason = seq + " - Failed reason: Read that did not appear in the correct orientation somewhere in the set of contig"
            seq_not_belong.append(failed_reason)
        elif contig_reads[seq]['contig_a'] is None or contig_reads[seq]['contig_b'] is None:
            failed_reason = seq + " - Failed reason: read from 1 of the contigs but do not span gaps between the contigs"
            seq_not_belong.append(failed_reason)
        else:
            if contig_reads[seq]['contig_a'] == contig_reads[seq]['contig_b']:
                distance = contig_reads[seq]['end_b'] - contig_reads[seq]['start_a']
                if distance < 1990 or distance > 2010:
                    failed_reason = seq + " - Failed reason: the distance from the beginnig of the first read to the end of the second read is not 2000 ± 10."
                    seq_not_belong.append(failed_reason)

    for each in seq_not_belong:
        print(each)

    '''print_str = ''
    for i in range(len(seq_not_belong) - 1):
        print_str += seq_not_belong[i] + ', '
    print_str += seq_not_belong[len(seq_not_belong) - 1]
    print(f"The following sequences have reads that don’t belong to the assembly: {print_str}")'''

    # -  The read pairs where both reads match to (distinct) contigs in supercontig1,
    #   and their estimated gap lengths
    seq_diff_contig = []
    for seq in contig_reads:
        if contig_reads[seq]['contig_a'] is not None and contig_reads[seq]['contig_b'] is not None:
            if contig_reads[seq]['contig_a'] != contig_reads[seq]['contig_b']:
                seq_diff_contig.append(seq)

    for sequence in seq_diff_contig:
        print(sequence + ': ' + str(contig_reads[sequence]))

    # report the number of mated pairs of reads
    # where the first read (the a read) maps to one contig and the second read (the b read) maps to another
    supercontig1 = []
    for sequence in seq_diff_contig:
        if contig_reads[sequence]['contig_a'] == 'contig_1':
            supercontig1.append(sequence)
    num_mated_pairs = len(supercontig1)
    print(f"The number of mated pairs of reads——where the first read (the a read) maps to one contig and the second read (the b read) maps to another——is {num_mated_pairs}")

    # report the possible length of the gap between the
    # two contigs (lower and upper bound)
    contig_1_length = len(contigs['contig1'])
    list = []
    for sequence in supercontig1:
        list_seq = []
        list_seq.append(sequence)
        first_read = contig_1_length - contig_reads[sequence]['start_a'] + 1
        second_read = contig_reads[sequence]['end_b']
        upper_bound = 2000 + 10 - first_read - second_read
        lower_bound = 2000 - 10 - first_read - second_read
        size_middle = (upper_bound + lower_bound)/2
        range_ = upper_bound - size_middle
        size_range = str(size_middle) + ' +/- ' + str(range_)
        # append the read of the values - since list can't append multiple values at once
        list_seq.append(size_range)
        list_seq.append(lower_bound)
        list_seq.append(upper_bound)
        list.append(list_seq)
    for item in list:
        print(str(item[0]) + ' ' + str(item[1]) + ' ' + str(item[2]) + ' ' + str(item[3]))


def find_contig_reads(reads, contigs):
    """
    Determine the contig in which each sequencing read appears (if any), along with
    where in the contig it matches.  The `contigs` dict will have a collection of 
    contig names mapped to contig sequences.  The `reads` dict contains separate keys 
    (labeled 'a' and 'b') for the two reads in each read pair.  However, this function
    returns a dictionary with a single key for each read pair (i.e., without 'a' or 'b').

    The value for each read-pair key will itself be a dictionary with the following keys:

        - 'contig_a' (str): the contig in which read 'a' was found, as 
          <contig name> or None (<contig name> will be a key in `contigs` dict)
        - 'start_a' (int): the start position (1-indexed) read 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'end_a' (int): the end position (1-indexed) read 'a' mapped to
          within its respective contig (None if not found in any contig)
        - 'contig_b' (str): the contig in which read 'b' was found, as
          <contig name> or None. (see 'contig_a' for example).
        - 'start_b' (int): the start position (1-indexed) read 'b' mapped to
          within its respective contig (None if not found in any contig)
        - 'end_b' (int): the end position (1-indexed) read 'b' mapped to
          within its respective contig (None if not found in any contig)

    The returned dictionary should look something like:
    {
        'seq1': {
            'contig_a': 'contig1',
            'start_a': 301,
            'end_a': 800,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        'seq2': {
            'contig_a': 'contig2',
            'start_a': 1101,
            'end_a': 1600,
            'contig_b': 'contig1'
            'start_b': 201,
            'end_b': 700,
        },
        'seq3' : {
            'contig_a': None,
            'start_a': None,
            'end_a': None,
            'contig_b': None
            'start_b': None,
            'end_b': None,
        },
        ...
    }

    Arguments:
        reads (dict str to str): dictionary of reads, mapping read names to sequences
        contigs (dict str to str): dictionary of contigs, mapping contig names to sequences

    Returns:
        Dictionary mapping read-pairs to information about their reads' locations in contigs.
    """

    #
    # YOUR CODE HERE
    #

    returned_dict = {}
    count = 1
    # loop over every read
    for read in reads:
        number = count
        found = False
        # compare every read to every contig
        for contig in contigs:
            if reads[read] in contigs[contig]:
                contig_number = contig[-1]
                contig_name = 'contig_' + str(contig_number)
                # report start and end
                start = contigs[contig].find(reads[read]) + 1
                end = start + len(reads[read])
                input_key = 'seq' + str(number)
                # capitalize on the fact that if a key is not already in the dict, then it means it's contig_a
                # if it's already in the dict, then it's contig_b
                if input_key not in returned_dict:
                    returned_dict[input_key] = {'contig_a': contig_name, 'start_a': start, 'end_a': end}
                    count = count
                else:
                    returned_dict[input_key]['contig_b'] = contig_name
                    returned_dict[input_key]['start_b'] = start
                    returned_dict[input_key]['end_b'] = end
                    count += 1
                found = True
        # if a read is matched with any of the contig, we don't need to check the rest of the contigs in the dict
        if not found:
            contig_name = None
            start = None
            end = None
            input_key = 'seq' + str(number)
            # similar logic when read-contig is matched
            if input_key not in returned_dict:
                returned_dict[input_key] = {'contig_a': contig_name, 'start_a': start, 'end_a': end}
                count = count
            else:
                returned_dict[input_key]['contig_b'] = contig_name
                returned_dict[input_key]['start_b'] = start
                returned_dict[input_key]['end_b'] = end
                count += 1

    return returned_dict


if __name__ == '__main__':
    """Call run_assembly_tester(). Do not modify this block."""
    run_assembly_tester()
