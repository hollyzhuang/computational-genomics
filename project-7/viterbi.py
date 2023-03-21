import sys
from math import log
from compsci260lib import get_fasta_dict


def run_viterbi():
    hmm_file = 'HMM.methanococcus.txt'
    input_file = 'bacterial.genome.fasta'
    print("Decoding the Viterbi path for %s using %s" % (input_file, hmm_file))

    vit_ret = viterbi_decoding(input_file, hmm_file)

    # Collect the segment counts for each state
    counts = count_segments(vit_ret)

    # Report the first 10 and last 10 segments of your decoding
    #
    # YOUR CODE HERE
    #
    for i in range(10):
        print(f"Segment {i+1}: {vit_ret[i]}")
    total_counts = counts['W'] + counts['S']
    for j in range(total_counts-10, total_counts):
        print(f"Segment {j+1}: {vit_ret[j]}")

    # Then report the number of segments that exist in each state.
    #
    # YOUR CODE HERE
    #
    print(f"The number of segments that exist in W state - AT-rich - is {counts['W']}.")
    print(f"The number of segments that exist in S state - GC-rich - is {counts['S']}.")


def viterbi_decoding(input_file, hmm_file):
    """Calculate the Viterbi decoding of an input sequence

    Arguments:
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
    
    print("Read a sequence of length", len(emit_str))
    
    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)]

    # Initialize the first column of the matrix
    for k in range(K):
        in_index = emitted_symbols[emit_str[0].upper()]
        viterbi[k][0] = emit_probs[k][in_index] + initial_probs[k]
    
    # Build the matrix column by column
    for i in range(1, len(emit_str)):
        in_index = emitted_symbols[emit_str[i].upper()]
        
        for l in range(K):
            # Compute the entries viterbi[l][i] and pointers[l][i]
            # Tip: Use float('-inf') for the value of negative infinity
            #
            # YOUR CODE HERE
            #
            max_store = []
            for j in range(K):
                # update the equation
                var = viterbi[j][i-1] + transitions[j][l] + emit_probs[l][in_index]
                # store the value in max_store
                max_store.append(var)
            # take the maximum value
            viterbi[l][i] = max(max_store)
            # new input in the dp table
            pointers[l][i] = max_store.index(max(max_store))

    # Traceback, stored as a list of segments in each state (represented using dictionaries)
    #
    # YOUR CODE HERE
    #
    return_path = [-1 for _ in range(len(emit_str))]

    # extract the maximum value of the last column
    last_column = []
    for i in range(K):
        last_column.append(viterbi[-1][i])
    last_index = last_column.index(max(last_column))
    return_path[-1] = last_index

    # start from the last column of the pointer matrix and end at the second column because the first column of the
    # pointer matrix did not have any input value (default zero)
    for j in reversed(range(0, len(emit_str) - 1)):
        return_path[j] = pointers[return_path[j+1]][j+1]

    # go through the indexes
    return_dict = []
    state = -1
    # convert into states and calculate start and end
    for k in range(len(return_path)):
        current_state = return_path[k]
        if current_state != state:
            start = k + 1
            state = current_state
        if k == len(return_path) - 1 or return_path[k] != return_path[k+1]:
            end = k + 1
            return_dict.append({'start': start, 'end': end, 'state': states[return_path[k]]})

    return return_dict


def count_segments(vit_ret):
    """Calculate the number of segments appearing in each state of
    the viterbi path

    Arguments:
        vit_ret (list of dicts): dictionary of segments in each state.
            see: return value of viterbi_decoding

    Returns:
        a dictionary mapping states to number of segments of that state. 
        e.g. {'W': 10, 'S': 9}
    """

    #
    # YOUR CODE HERE
    #

    # initialize the values
    new_dict = {}
    new_dict['W'] = {}
    new_dict['S'] = {}
    new_dict['W']['count'] = 0
    new_dict['S']['count'] = 0
    new_dict['W']['length_dist'] = {}
    new_dict['S']['length_dist'] = {}

    # the following code will be able to output the number of states (W and S) as well as the length distribution
    # of the segments of each state
    # the final output only prints out the number od states to pass test cases in gradescope
    for seg in vit_ret:
        if seg['state'] == 'W':
            new_dict['W']['count'] += 1
            length = seg['end'] - seg['start'] + 1
            if length not in new_dict['W']['length_dist']:
                new_dict['W']['length_dist'][length] = 1
            else:
                new_dict['W']['length_dist'][length] += 1
        if seg['state'] == 'S':
            new_dict['S']['count'] += 1
            length = seg['end'] - seg['start'] + 1
            if length not in new_dict['S']['length_dist']:
                new_dict['S']['length_dist'][length] = 1
            else:
                new_dict['S']['length_dist'][length] += 1

    return_dict = {'W': new_dict['W']['count'], 'S': new_dict['S']['count']}
    return return_dict


if __name__ == '__main__':
    """Call run_viterbi(), do not modify"""
    run_viterbi()
