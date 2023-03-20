import sys
import random
from compsci260lib import *

def run_simulate():
    """
    Simulates the sequencing process then empirically compute and report 
    the quantities from parts a-c.
    """

    iterations = 20
    G = 3000000
    R = 40000
    L = 400

    # Call simulate(G, R, L) `iteration` times. Report the empirical values 
    # from parts a-c for each iteration
    #
    # YOUR CODE HERE
    #
    empirical_iteration_output = []
    # for every iteration, run the simulation
    for i in range(iterations):
        tuple_output = simulate(G, R, L)
        empirical_iteration_output.append(tuple_output)

    print(empirical_iteration_output)

    # report the output in a readable manner
    for i in range(len(empirical_iteration_output)):
        iteration = i+1
        interation_output = empirical_iteration_output[i]

        print(f"For iteration {iteration}, the coverage is {interation_output[0]}, the number of nucleotides not covered is {interation_output[1]}, the number of contigs is {interation_output[2]}, and the average length of contigs is {interation_output[3]}.")

    # take the average value for each value in the tuple over the 20 iterations
    # store the sum of each variable so that they can be divided by the iteration numbers to calculate average
    sum_coverage = 0
    sum_uncovered_nucleotides = 0
    sum_contigs = 0
    sum_contigs_len = 0
    # iterate through every iteration and add up each individual value corresponding to each other from tuples
    for iteration in empirical_iteration_output:
        sum_coverage += iteration[0]
        sum_uncovered_nucleotides += iteration[1]
        sum_contigs += iteration[2]
        sum_contigs_len += iteration[3]
    # take the average for each value
    ave_coverage = sum_coverage/iterations
    ave_uncovered_nucleotides = sum_uncovered_nucleotides/iterations
    ave_contigs = sum_contigs/iterations
    ave_contigs_len = sum_contigs_len/iterations

    print(ave_coverage, ave_uncovered_nucleotides, ave_contigs, ave_contigs_len)

def simulate(G, R, L):
    """
    Simulates one iteration of the sequencing process and empirically compute 
    the empirical coverage (average number of times a nucleotide in the genome
    was sequenced), the number of nucleotides not covered by any read, the 
    number of contigs assuming you can use the oracular assembly algorithm to
    assemble all the reads, and the average length of these contigs.

    Args:
        G (int) - the length of the genome
        R (int) - the number of reads
        L (int) - the length of each read

    Returns
        a tuple of floats:

            (Empirical coverage, 
             Number of nucleotides not covered by any read, 
             Number of contigs,
             Average length of these contigs)
    """
    #
    # YOUR CODE HERE
    #
    # initialize all the elements in the list of length G to 0
    list_count = [0] * G
    # iterate the process for R number of reads
    for i in range(R):
        # randomly select a starting point
        start = random.randrange(0, G)
        # based on the starting point, locate the ending point because read length is the same across
        end = start + L
        # the following code accounts for the edge effect, so i solve it by connecting the end of the sequence to the
        # start of the sequence. whenever end is larger than the total length of G, the read match extends to the front
        # of the genome
        if end <= G:
            for j in range(start, end):
                list_count[j] += 1
        # when end is larger than G
        else:
            exceed = end - G
            within = L - exceed
            # first half of the read match
            for k in range(start, start + within):
                list_count[k] += 1
            # second half of the read match that starts from the front of the genome
            for t in range(exceed):
                list_count[t] += 1

    # calculate coverage: add every number up in every location of the genome and divided by the length of the genome
    sum_nucleotides = 0
    for nucleotide in list_count:
        sum_nucleotides += nucleotide
    coverage = sum_nucleotides / G

    # calculate the number of nucleotides not covered by any read: number of zeros in the list
    zero_count = 0
    for nucleotide in list_count:
        if nucleotide == 0:
            zero_count += 1

    # calculate the number of contigs: whenever there's zero/zeros behind non-zero numbers, count as 1 contig
    contig_count = 0
    i = 0
    last_not_zero = False
    while i < G:
        # if the index position in the genome is not zero, last_not_zero is true, and proceed to the next position
        if list_count[i] != 0:
            last_not_zero = True
            # proceed to the next position
            i += 1
        # when the index position in the genome is zero and if last_not_zero is true, count as 1 contig
        else:
            if last_not_zero:
                contig_count += 1
                # when found, last_not_zero become false
                last_not_zero = False
                # proceed to the next position
                i += 1
            # the case of consecutive zeros: when the next zero appears after zeros, it's noa contig
            else:
                # proceed to the next position
                i += 1

    #  calculate average length of contigs: similar logic to calculating the number of contigs
    store_length = []
    i = 0
    last_not_zero = False
    first_index = 0
    while i < G:
        if list_count[i] != 0:
            # first non-zero appearance
            if not last_not_zero:
                first_index = i
            # not the first non-zero appearance, proceed
            last_not_zero = True
            i += 1
        else:
            # first appearance of zero after series of non-zeros
            if last_not_zero:
                last_index = i
                # get the length by last - first
                contig_length = last_index - first_index
                store_length.append(contig_length)
                # last not zero is no long true
                last_not_zero = False
                i += 1
            # consecutive zeros
            else:
                i += 1

    # get the sum of all the contig length
    sum_length = 0
    for length in store_length:
        sum_length += length

    # get the average length by dividing the sum length with contig numbers
    aveg_length = sum_length / contig_count

    return (coverage, zero_count, contig_count, round(aveg_length, 4))

if __name__ == '__main__':
    """Call run_simulate(), do not modify"""
    run_simulate()
