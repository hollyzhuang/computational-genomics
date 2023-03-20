from bwt import *
from fm_index import *
from compsci260lib import *


def align_patient_reads():
    """Align patient reads to each bacterial genome"""

    # Patients
    patients = ['patient1',
                'patient2',
                'patient3']

    # Bacterial species
    panel = ['Bacteroides_ovatus',
             'Bacteroides_thetaiotaomicron',
             'Bifidobacterium_longum',
             'Eubacterium_rectale',
             'Lactobacillus_acidophilus',
             'Peptoniphilus_timonensis',
             'Prevotella_copri',
             'Roseburia_intestinalis',
             'Ruminococcus_bromii',
             'Vibrio_cholerae']

    # Store the genome sequence and bwt structures for each bacterial species
    bact_sequences = {}
    bact_bwt_structures = {}
    for bacteria in panel:
        # For each of the bacteria in our panel, create a dictionary entry
        # where keys are bacterial names and values are their genome sequences.
        bact_sequences[bacteria] = list(get_fasta_dict(
            'reference_genomes/%s.fasta' % bacteria).values())[0]

        # For each of the bacterial genome sequences, create a dictionary entry
        # where keys are bacterial names and values are their BWT data
        # structures.
        bact_bwt_structures[bacteria] = make_all(bact_sequences[bacteria])

    # Store a mapping of patient names to reads
    patient_reads = {}
    for patient in patients:
        reads = list(get_fasta_dict('patients/%s.fasta' % patient).values())
        patient_reads[patient] = reads[:100]

    # print(patient_reads)
    # Store the reads mapped per bacterial species per patient
    mapped_patient_reads = {}

    # Consider all patients
    for patient in patients:
 
        # Find uniquely mapped reads for each bacteria for the patient
        # and store the read start positions
        mapped_reads = find_aligned_reads_for_patient(
            patient_reads[patient], bact_bwt_structures)
        mapped_patient_reads[patient] = mapped_reads
    # print(mapped_patient_reads)

    # Report the microbe prevalences for each of your patients
    # 
    # YOUR CODE GOES HERE
    #
    # dictionary to store total mapped reads for each patient
    total_mapped_reads = {}
    # loop through every patient
    for patient in mapped_patient_reads:
        total_count = 0
        # after looping through every bacterium, calculate the sum of the reads that get matched while satisfying
        # all requirements to be considered as a valid read match
        for bacteria in mapped_patient_reads[patient]:
            total_count += len(mapped_patient_reads[patient][bacteria])
        # each patient has a total sum of the reads matched for 10 microbes
        total_mapped_reads[patient] = total_count
    print(total_mapped_reads)

    # look at the prevalence for EACH patient
    prevalence_dict = {}
    for patient in mapped_patient_reads:
        patient_total = total_mapped_reads[patient]
        # a dictionary of each key's value as another dictionary
        prevalence_dict[patient] = {}
        for bacteria in mapped_patient_reads[patient]:
            # count number of matches for each microbe
            count = len(mapped_patient_reads[patient][bacteria])
            # prevalence = count number of matches for 1 microbe / total matches for 10 microbes
            prevalence = round(count/patient_total, 4)
            # store the prevalence value in the new dictionary
            prevalence_dict[patient][bacteria] = prevalence

    print(prevalence_dict)

    # Use `read_mapper` and `longest_zeroes` to identify unmapped regions
    # for the relevant patients and species (questions 2f-h)
    # 
    # YOUR CODE GOES HERE
    #
    # patient 1 and patient 2's mapped reads to bacteria 'Vibrio_cholerae'
    patient_1_mapped_reads = mapped_patient_reads['patient1']['Vibrio_cholerae']
    patient_2_mapped_reads = mapped_patient_reads['patient2']['Vibrio_cholerae']
    # store Vibrio_cholerae's genome
    vibrio_genome_seq = bact_sequences['Vibrio_cholerae']

    # for patient 1 (asymptomatic)
    count_vector_1 = read_mapper(patient_1_mapped_reads, vibrio_genome_seq)
    patient_1_return = longest_zeros(count_vector_1)
    print(patient_1_return)

    # for patient 2 (symptomatic)
    count_vector_2 = read_mapper(patient_2_mapped_reads, vibrio_genome_seq)
    patient_2_return = longest_zeros(count_vector_2)
    print(patient_2_return)

    # extract the corresponding sequence from the reference genome mapped with reads from patient 1
    extract_seq = vibrio_genome_seq[patient_1_return[0]:patient_1_return[1]+1]
    print(extract_seq)

def find_aligned_reads_for_patient(reads, bact_bwt_structures):
    """
    Given a list of reads for a patient, identify the reads that uniquely map
    to each bacteria's genome using its relevant BWT data structure. Reads that
    are mapped to a bacteria's genome should be stored as starting positions in
    a list.

    Args:
        reads (list of str):
            mapping a patient's read names (str) to read sequences (str).

        bact_bwt_structures: dictionary mapping bacterial names (str) to 
            structures required for efficient exact string matching). Refer
            to the return tuple from `make_all` in fm_index.py.

    Returns:
        a dictionary mapping bacterial names (str) to the start positions of
        reads mapped to that bacteria's genome. Note that the end positions are
        not needed as the reads are all 50 bases long. Start positions should
        be stored using 0-indexing.

        Example return structure:
        {
            'Bacteroides_ovatus': [8, 124, 179, ...],
            ...
        }
    """

    # return counts: for each patient
    counts = {}
    panel = list(bact_bwt_structures.keys())
    # 
    # YOUR CODE GOES HERE
    #
    # initialize the key values as names of bacteria in the return dictionary counts
    for bacteria in panel:
        counts[bacteria] = []

    # loop through every read
    for read in reads:
        # Reverse complement the read first and then match it to the genome is more simple and efficient: 1). index tracking is
        # not going to change with the reverse complement of the read; 2). the read is way smaller than the genome, so reverse
        # complement the read is much faster than reverse complementing the whole genome.
        read = reverse_complement(read, seq_type='DNA')
        # keep track how many times a read is matched to the 10 microbes' genome; if more than 1, it means one read is
        # matched to more than one genome, and we don't consider that read when calculating prevalence
        match_genome_count = 0
        # store the bacterium's name of the genome that's matched
        unique_bacterial_genome = ""
        # store index
        index_store = 0
        for bacteria in panel:
            bacteria_sequence_bwt_data = bact_bwt_structures[bacteria]
            # call the function find() on each read for each bacterial genome
            match_index_list = find(read, bacteria_sequence_bwt_data)
            # check if match_index_list is a list or a string
            # if it's a string, from the fm.index file, it means that no match is found as it will output "not found"
            # if it's a list, proceed, because at least one match is found
            if isinstance(match_index_list, list):
                # whenever a read match is found in any of the 10 bacterial genomes, match_genome_count increases by 1
                match_genome_count += 1
                unique_bacterial_genome = bacteria
                # when multiple match sequences are found in a genome for a read, only record the first one
                # keep it as evidence of that bacterium's prevalence, but count it only once
                # we want to check if a read appears in a bacterial genome, so recording one location is enough to show
                # that it's present, and for calculating prevalence it also makes more sense - it doesn't matter how
                # many times the read can get match in a single bacterial genome since at least one site can match
                # to the read
                index_store = match_index_list[0]
        # add to the dictionary only when one match is found (not more than 1 or 0)
        if match_genome_count == 1:
            counts[unique_bacterial_genome].append(index_store)
    # sort the index
    for count in counts:
        counts[count] = sorted(counts[count])
    return counts


def read_mapper(starting_positions, genome):
    """
    Using the starting positions of reads that were aligned to a bacterial
    genome, construct a count vector (as a list) that counts how many reads
    were aligned to a position in the genome. The vector's size will be the 
    length of the genome. You may assume that each read is 50 base pairs long.

    Args:
        starting_positions (list of ints): 
            starting positions of reads that were aligned to a genome.

        genome (str): the genomic nucleotide sequence to which the reads 
                      were aligned.

    Returns:
        (list of ints): vector of aligned read counts to the genome.
        i.e. [c_1, c_2, ..., c_i, ..., c_n], where n=length of the genome
        and c_i = the count of aligned reads for the patient at genome
        position i.
    """

    # 
    # YOUR CODE GOES HERE
    #
    length = len(genome)
    # create a count vector with length of the input genome
    count_vector = [0]*length
    # loop through every start index location
    for start in starting_positions:
        # each read is assumed to be the same length, 50 bases
        end = start + 50
        # increment by 1 for the read match corresponding to the index of the count vector
        for i in range(start, end):
            count_vector[i] += 1
    return count_vector

# It may be helpful to read the documentation for the methods
# given below, but you will NOT have to make any changes to
# them in order to complete the problem set.


def longest_zeros(count_vector):
    """Given a count vector, return the start and stop position (inclusive) of
    the longest string of internal zeros in the vector. If there is no
    internal string of zeros, return None.

    Examples to understand behavior:
    input -> output
    contain valid internal string of zeros:
    [1, 1, 1, 0, 0, 1, 1] -> (3, 4)
    [1, 1, 1, 0, 1, 1, 1]  -> (3, 3)
    [0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0] -> (5, 6)

    do not contain internal string of zeros:
    [0, 0, 0, 1, 1, 1, 1] -> None
    [1, 1, 1, 1, 0, 0, 0] -> None
    [0, 0, 0, 0, 0, 0, 0] -> None
    [1, 1, 1, 1, 1, 1, 1] -> None

    Args:
        count_vector (list of ints): vector of aligned read counts to the
        genome.  see: the return value for `read_mapper()`

    Returns:
        (tuple of (int, int)): the start and stop position (inclusive) of the longest
        internal string of zeros in the count_vector. If there are no internal
        runs-of-zero the return will be None.
    """

    zero_nums = []
    genome_len = len(count_vector)
    for i in range(0, genome_len):
        if count_vector[i] == 0:
            zero_nums.append(i)

    if len(zero_nums) == 0:
        return None

    counter = 1
    longest_run = {}
    longest_run[1] = []
    longest_run[1].append(zero_nums[0])
    for z in zero_nums[1:]:
        if (z - longest_run[counter][-1]) == 1:
            longest_run[counter].append(z)
        else:
            counter += 1
            longest_run[counter] = []
            longest_run[counter].append(z)

    for run in list(longest_run.keys()):
        if 0 in longest_run[run] or genome_len-1 in longest_run[run]:
            del longest_run[run]

    longest = []
    for run in list(longest_run.values()):
        if len(run) > len(longest):
            longest = run

    # There was no run of zeroes found, return None
    if len(longest) == 0:
        return None

    start = longest[0]
    stop = longest[-1]

    # Return the start and end positions of the longest zeroes (0-indexed)
    return start, stop


if __name__ == '__main__':
    align_patient_reads()
