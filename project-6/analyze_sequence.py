from generate_HMM_sequence import generate_HMM_sequence


def run_analyze_sequence(gen_seq_file):
    """Load the nucleotide sequence HMM file and generate a 100000 length sequence,
    then do some analysis of the results, including the average length in each state.

    Args:
        gen_seq_file: The name of the output file for the generated nucleotide sequence.
    """
    seq_length = 100000
    state_sequence, observed_sequence = generate_HMM_sequence("HMM.parameters.txt", seq_length)

    # Write your generated sequence to file. Report the states that prevailed whenever
    # TCGA was observed, and their frequencies when observing TCGA. Finally, report
    # the average length spent in each state across the entire sequence.
    #
    # YOUR CODE GOES HERE
    #

    # write the generated sequence to file
    state_str = ''
    seq_str = ''
    # convert from lists to two strings for display in the output file
    for state in state_sequence:
        state_str = state_str + state
    for seq in observed_sequence:
        seq_str = seq_str + seq
    filename = "nucleotide_sequence.txt"
    # write the file
    f = open(filename, "w")
    # two strings seperated by a line
    f.write(state_str + " \n" + seq_str)
    f.close()

    # the states that prevailed whenever
    # TCGA was observed, and their frequencies when observing TCGA
    subsequence_list = ['T', 'C', 'G', 'A']
    freq_dict = compute_state_frequencies(state_sequence, observed_sequence, subsequence_list)
    print(freq_dict)
    for state in freq_dict:
        print(f"When the observed nucleotide sequence was TCGA, the frequency of the hidden state sequence {state} is {freq_dict[state]}.")

    # the average length spent in each state across the entire sequence
    length_dict = compute_average_length(state_sequence)
    print(length_dict)
    for length in length_dict:
        print(f"The average length for all-{length} state sequence is {length_dict[length]}")


def compute_average_length(state_sequence):
    """Given a state sequence, return the average length in each state

    Arguments:
        state_sequence: generated sequence of states, represented as list of single-char strings

    Returns:
        a dictionary mapping the state name to the average length in the state.
        Example return:
        {
            "W": 5.5,
            "S": 4.5
        }
    """
    # Check the type of state_sequence
    if type(state_sequence) is not list:
        raise TypeError(f"The argument 'state_sequence' must be a list of strings. "
                        f"(received type {type(state_sequence)})")

    # Compute the average length in each state
    #
    # YOUR CODE HERE
    #
    # initialize the count and total lengths for W and S
    total_length_w = 0
    total_length_s = 0
    count_w = 0
    count_s = 0
    state = ''

    # loop through the whole state sequence to count for W and S
    for i in range(len(state_sequence)):
        # record the current state
        current_state = state_sequence[i]
        # if current state is not the previous state that has remained unchanged for a while
        # record the start location and change the state
        if current_state != state:
            start = i
            state = current_state
        # if state_sequence[i] is the last in the state sequence list or the next in
        # sequence does not equal to the current state (meaning that it's the end of a sequence of continuous W or S)
        # record end position and calculate the length
        if i == len(state_sequence) - 1 or state_sequence[i] != state_sequence[i+1]:
            end = i
            length = end - start + 1
            # based on whether it's W or S
            # record their count and add up the total length
            if state_sequence[i] == 'W':
                count_w += 1
                total_length_w += length
            if state_sequence[i] == 'S':
                count_s += 1
                total_length_s += length

    # create a dictionary to store value and output them
    dict_ave_length = {'W': total_length_w / count_w, 'S': total_length_s / count_s}

    return dict_ave_length


def compute_state_frequencies(state_sequence, observed_sequence, subsequence):
    """Given the state and observed sequences, return the state sequences that
    emitted the query subsequence and frequency in which those state sequences
    sequences emitted the subsequence.

    Arguments:
        state_sequence (list of single-char strings): generated state sequence
        observed_sequence (list of single-char strings): generated observed sequence
        subsequence (list of single-char strings): the observed subsequence to count

    Returns:
        a dictionary mapping the state name to the frequency of observing the
        provided sequence. Example return:
        {
            "WWWW": 2,
            "WWWS": 1,
            ...
        }
    """

    # Check the types for each of the input arguments
    if type(state_sequence) is not list:
        raise TypeError(f"The argument 'state_sequence' must be a list of strings. "
                        "(received type {type(state_sequence)})")
    if type(observed_sequence) is not list:
        raise TypeError(f"The argument 'observed_sequence' must be a list of strings. "
                    "(received type {type(observed_sequence)})")
    if type(subsequence) is not list:
        raise TypeError(f"The argument 'subsequence' must be a list of strings. "
                    "(received type {type(subsequence)})")

    #
    # YOUR CODE HERE
    #

    # initialize the dictionary
    return_dict = {}
    # check if match
    match = False
    # make code reproducible no matter what the length of the subsequence is
    sub_length = len(subsequence)
    # loop over every group of subsequence on the observed_sequence list
    for i in range(len(observed_sequence)-sub_length+1):
        previous_match = True
        for j in range(sub_length):
            # only if previous match is true - the program proceed
            if previous_match:
                # if characters match, then the previous_match condition stays true
                if observed_sequence[i+j] == subsequence[j]:
                    match = True
                    previous_match = match
                # if characters don't match, go back to the for loop for i
                else:
                    match = False
                    previous_match = match
        # if all matches with the subsequence, then find the corresponding state and put into the dictionary
        if match:
            state_seq = ''
            # get the corresponding state
            for k in range(sub_length):
                state_seq += state_sequence[i+k]
            # check if the found state is already in the dictionary
            if state_seq in return_dict:
                return_dict[state_seq] += 1
            else:
                return_dict[state_seq] = 1
    return return_dict


if __name__ == "__main__":
    """Main method call, do not modify"""
    run_analyze_sequence("nucleotide_sequence.txt")
