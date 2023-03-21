import os
import random


def generate_HMM_sequence(hmm_file, seq_length):
    """Load an HMM specification from file, and generate state and observed
    sequences through sampling of the HMM.

    The HMM will have states labeled as strings and observations as characters
    which will be useful when we generate an HMM with nucleotide sequence
    observations.

    An example return sequence may look like:
        (["W", "S", "W"], ["A", "T", "G"]) or
        (["F", "L", "L"], ["2", "6", "6"])

    Arguments:
        hmm_file (str): path to the HMM file
        seq_length (int): the length of the sequence we will generate

    Returns:
        a tuple containing two lists of strings, with the first list storing
        the sequence of states (which are arbitrary strings) and the second
        list storing the sequence of observations (which are assumed to be
        single-character strings):

        (state sequence as a list of strings,
         observed sequence as a list of single-character strings)
    """

    if not os.path.exists(hmm_file):
        print("Can't open HMM parameter file: %s" % hmm_file)
        return -1

    f = open(hmm_file, "r")

    # read the state names
    states = f.readline().strip().split()

    # read the initial probabilities
    initial_probs = f.readline().strip().split()
    initial_probs = [float(p) for p in initial_probs]

    # read the transition matrix
    transitions = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        transitions[state] = [float(p) for p in matrix_row]

    # read the input alphabet
    input_alphabet = f.readline().strip().split()

    # read the emission matrix
    emission = {}
    for i in range(0, len(states)):
        state = states[i]
        matrix_row = f.readline().strip().split()
        emission[state] = [float(p) for p in matrix_row]
        # normalize
        sum_emission = sum(emission[state])
        emission[state] = [x / sum_emission for x in emission[state]]
    f.close()

    #
    # YOUR CODE HERE
    #

    # initialize the two lists that will be returned: state sequence list and dna list (ATCG...)
    state_seq = []
    dna_seq = []

    # initialize the situation for the initial probabilities for W and S
    # get the first state based on initial_probs
    first_state = random.choices(states, weights=initial_probs, k=1)
    first_state = first_state[0]

    # get the first sequence based on whether the first state is W or S
    first_seq = []
    if first_state == 'W':
        first_seq = random.choices(input_alphabet, weights=emission['W'], k=1)
    if first_state == 'S':
        first_seq = random.choices(input_alphabet, weights=emission['S'], k=1)
    first_seq = first_seq[0]

    # append to both of the lists
    state_seq.append(first_state)
    dna_seq.append(first_seq)

    # set previous to first
    previous_state = first_state

    # generate sequence for seq_length times
    for _ in range(seq_length):
        # determine the state
        next_state = []
        if previous_state == 'W':
            next_state = random.choices(states, weights=transitions['W'], k=1)
        if previous_state == 'S':
            next_state = random.choices(states, weights=transitions['S'], k=1)
        next_state = next_state[0]

        # determine which dna sequence based on the state
        next_seq = []
        if next_state == 'W':
            next_seq = random.choices(input_alphabet, weights=emission['W'], k=1)
        if next_state == 'S':
            next_seq = random.choices(input_alphabet, weights=emission['S'], k=1)
        next_seq = next_seq[0]

        # append them to both lists
        state_seq.append(next_state)
        dna_seq.append(next_seq)

        # set previous state to next
        previous_state = next_state

    return state_seq, dna_seq  # a tuple containing the state and observation sequences
