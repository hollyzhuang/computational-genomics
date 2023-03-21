import sys
from math import log, exp
from compsci260lib import get_fasta_dict


def run_posterior():
    input_file = "bacterial.genome.fasta"
    hmm_file = "HMM.methanococcus.txt"

    posterior = posterior_decoding(input_file, hmm_file)

    # Report the first and last ten segments in your decoding
    #
    # YOUR CODE HERE
    #
    for i in range(10):
        print(f"Segment {i+1}: {posterior[i]}")
    for j in range(len(posterior)-10, len(posterior)):
        print(f"Segment {j+1}: {posterior[j]}")

    # to check what tRNAs are presented
    print(posterior[50:101])


def posterior_decoding(input_file, hmm_file):
    """
    Calculate the posterior decoding and return the decoded segments.

    input_file (str): path to input fasta file
    hmm_file (str): path to HMM file

    Returns:
        A list of dictionaries of segments in each state (1-indexed). 
        An example output may look like:

        [
            {‘start’: 1, ‘end’: 12, ‘state’: ‘S’},
            {‘start’: 13, ‘end’: 20, ‘state’: ‘W’},
            ...
        ]

    """

    # Open HMM description file
    try:
        f_hmm_file = open(hmm_file, 'r')
    except IOError:
        print("IOError: Unable to open HMM file: %s." % hmm_file)
        print("Exiting.")
        sys.exit(1)

    # Read the state names
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # Read the initial probability vector (and log transform entries)
    probs = f_hmm_file.readline().split()
    initial_probs = [log(float(prob)) for prob in probs]
        
    # Read the transition matrix (and log transform entries)
    transitions = [None for _ in range(K)]
    for k in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(trans_prob)) for trans_prob in matrix_row_arry]
        transitions[k] = matrix_row
        
    # Read the emitted symbols
    emitted_symbols_list = f_hmm_file.readline().split()
    emitted_symbols = {symbol: index for index, symbol in enumerate(emitted_symbols_list)}
    
    # Read the emission probability matrix (and log transform entries)
    emit_probs = [None for _ in range(K)]
    for k in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [log(float(emit_prob)) for emit_prob in matrix_row_arry]
        emit_probs[k] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = list(seq_dict.values())[0]  # there's only 1
    
    print("Done reading sequence of length ", len(emit_str))
    
    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emitted_symbols,
        emit_probs, emit_str)

    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions, 
        emitted_symbols, emit_probs, emit_str)

    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[0 for _ in range(len(emit_str))] for _ in range(K)]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[0][i] versus posterior[1][i]
        for k in range(K):
            posterior[k][i] = forward[k][i] + backward[k][i]
    
    # Create the list of decoded segments to return
    #
    # YOUR CODE HERE
    #
    # create a sequence of most likely states in a list
    path_list = []
    for i in range(len(emit_str)):
        temporary_list = []
        for k in range(K):
            temporary_list.append(posterior[k][i])
        path_list.append(temporary_list.index(max(temporary_list)))

    # create the segment - similar code to viterbi
    return_list = []
    state = -1
    for k in range(len(path_list)):
        current_state = path_list[k]
        if current_state != state:
            start = k + 1
            state = current_state
        if k == len(path_list) - 1 or path_list[k] != path_list[k + 1]:
            end = k + 1
            return_list.append({'start': start, 'end': end, 'state': states[path_list[k]]})

    return return_list


def calc_log_sum(log_list):

    # order the values from large to small
    log_list.sort(reverse=True)

    # constant term
    log_sum = log_list[0]
    log_term = 1

    # apply the equation that can avoid underflow
    for i in range(1, len(log_list)):
        log_term += exp(log_list[i]-log_list[0])
    log_sum += log(log_term)

    return log_sum


def run_forward(states, initial_probs, transitions, emitted_symbols, 
                emit_probs, emit_str):
    """Calculates the forward (log) probability matrix.

    Arguments:
        states (list of str): list of states as strings
        initial_probs (list of float): list of log(initial probabilities) for each
            state
        transitions (list of list of float): matrix of log(transition probabilities)
        emitted_symbols (dict {str: int}): dictionary mapping emitted symbols to their index
        emit_probs (list of list of float): matrix of log(emission probabilities)
            for each state and emission symbol
        emit_str (str):

    Returns:
        (list of list of floats): matrix of forward (log) probabilities
    """

    K = len(states)
    N = len(emit_str)
    
    forward = [[0 for _ in range(N)] for _ in range(K)]
    
    # Initialize - same initialization for viterbi decoding
    emit_index = emitted_symbols[emit_str[0].upper()]
    for k in range(K):
        forward[k][0] = initial_probs[k] + emit_probs[k][emit_index]

    # Iterate
    for i in range(1, N):
        emit_index = emitted_symbols[emit_str[i].upper()]

        # Compute the forward probabilities for the states
        #
        # YOUR CODE HERE
        #
        for j in range(K):
            sum_store = []
            for k in range(K):
                # update the equation
                var = forward[k][i-1] + transitions[k][j] + emit_probs[j][emit_index]
                sum_store.append(var)
            # append the value to the dp table
            forward[j][i] = calc_log_sum(sum_store)

    return forward


def run_backward(states, initial_probs, transitions, emitted_symbols,
                 emit_probs, emit_str):
    """Calculates the backward (log) probability matrix.

        Arguments:
            states (list of str): list of states as strings
            initial_probs (list of float): list of log(initial probabilities) for
                each state
            transitions (list of list of float): matrix of log(transition
                probabilities)
            emitted_symbols (dict {str: int}): dictionary mapping emitted symbols to their index
            emit_probs (list of list of float): matrix of log(emission
                probabilities) for each state and emission symbol
            emit_str (str):

        Returns:
            (list of list of floats): matrix of backward (log) probabilities
    """

    
    K = len(states)
    N = len(emit_str)

    backward = [[0 for _ in range(N)] for _ in range(K)]

    # Initialize
    for k in range(K):
        backward[k][N - 1] = log(1)  # which is zero, but just to be explicit...

    # Iterate
    for i in range(N - 2, -1, -1):
        emit_index = emitted_symbols[emit_str[i + 1].upper()]

        # Compute the backward probabilities for the states
        #
        # YOUR CODE HERE
        #
        for j in range(K):
            sum_store = []
            for k in range(K):
                # update the equation
                var = backward[k][i+1] + transitions[j][k] + emit_probs[k][emit_index]
                sum_store.append(var)
            # append value to the dp table
            backward[j][i] = calc_log_sum(sum_store)
    return backward


if __name__ == '__main__':
    """Call run_posterior(), do not modify"""
    run_posterior()
