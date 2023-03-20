import random
import timeit

# random_list needs to be available as a global variable
# for the timeit function to work properly
random_list = None


def brute_force(score_list):
    """Get the maximum similarity score using brute force.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """
    
    #
    # Your code goes here
    #
    # initiate the highest sum
    highest_sum = -100
    # every substring is the prefix of some suffix
    # loop through the whole list from offset 1 to offset n, basically loop through the suffix
    for i in range(len(score_list)):
        sum = 0
        # loop through every prefix of the suffix
        for j in range(len(score_list) - i):
            # add the new element to the total sum
            sum += score_list[i+j]
            # if the sum is larger, than update the highest sum value
            if sum > highest_sum:
                highest_sum = sum
    return highest_sum


def divide_conquer(score_list, left=None, right=None):
    """Get the maximum similarity score using divide and conquer.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """
    #
    # Your code goes here
    #
    # base case
    if left is None and right is None:
        left = 0
        right = len(score_list) - 1

    # base case
    if right == left:
        return score_list[left]

    # get the middle index
    mid = int((left + right)/2)

    # maximum sublist sum from the left
    # middle index is included
    max_left = -10000
    sum = 0
    for i in range(left, mid + 2):
        sum += score_list[i]
        if sum > max_left:
            max_left = sum

    # maximum sublist sum from the right,
    # middle index is not included
    max_right = -10000
    sum = 0
    for i in range(mid + 1, right + 1):
        sum += score_list[i]
        if sum > max_right:
            max_right = sum

    # get max total by comparing left and right sub-lists
    max_total = max(divide_conquer(score_list, left, mid), divide_conquer(score_list, mid + 1, right))
    # max_left + max_right for including the sub-array that includes the middle element
    return max(max_total, max_left + max_right)

def linear(score_list):
    """Get the maximum similarity score in linear time.

    Args: score_list (list): list of integer similarity scores

    Returns: (int) of the maximal computed score
    """

    #
    # Your code goes here
    #
    # initiate MAX_SO_FAR and MAX_INCLUDING_HERE
    MAX_SO_FAR = -10000
    MAX_INCLUDING_HERE = 0
    # loop through every element in the list
    for i in range(len(score_list)):
        # update the maximum sum by adding the current index i to the former max at index i-1
        MAX_INCLUDING_HERE = MAX_INCLUDING_HERE + score_list[i]
        # update MAX_INCLUDING_HERE
        MAX_INCLUDING_HERE = max(MAX_INCLUDING_HERE, score_list[i])
        # update MAX_SO_FAR if the new MAX_INCLUDING_HERE is larger
        MAX_SO_FAR = max(MAX_SO_FAR, MAX_INCLUDING_HERE)
    return MAX_SO_FAR


def run_collective_similarity():
    """Run collective similarity and collect timing for each algorithm.

    Note: there is no reporting in this problem, `brute_force`,
    `divide_conquer`, and `linear` will be evaluated for correctness.
    """

    # Declare random_list as the same global variable defined
    # on line 7 above for use in the timeit library below.
    global random_list

    # You can use this to test the correctness of your code by using
    # sample_list as an input to each function.  You could also consider
    # creating other sample lists as tests, including a test case where
    # the right answer of zero occurs when an empty list would be optimal.

    sample_list = [2, -3, -4, 4, 8, -2, -1, 1, 10, -5]

    brute_force(sample_list)
    divide_conquer(sample_list)
    print(linear(sample_list))

    # This part below is used to test the runtime of your code, an example is
    # given below for brute force algorithm with a random list of length 100.
    # You will have to measure the runtime of each algorithm on every input size
    # given in the problem set.

    # """
    allowed_scores = [i for i in range(-10, 11)]
    random_list = [random.choice(allowed_scores) for _ in range(100000000)]
    linear_runtime = timeit.timeit("linear(random_list)",
                                       setup="from __main__ import linear, random_list",
                                       number=1)
    print(linear_runtime)
    # """


if __name__ == "__main__":
    """Run run_collective_similarity(). Do not modify this code"""
    run_collective_similarity()
