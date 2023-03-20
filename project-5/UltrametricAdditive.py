from compsci260lib import *


def solve_ultrametric_additive():
    
    # Distance metrics for table 1 and table 2
    dist_1 = {"1,2" : 0.3, "1,3" : 0.7, "1,4" : 0.7,
              "2,3" : 0.6, "2,4" : 0.6,
              "3,4" : 0.6}

    dist_2 = {"1,2" : 0.9, "1,3" : 0.4, "1,4" : 0.6, "1,5" : 0.9,
              "2,3" : 0.9, "2,4" : 0.9, "2,5" : 0.4,
              "3,4" : 0.6, "3,5" : 0.9,
              "4,5" : 0.9}

    # Check if dist_1 and dist_2 are ultrametric and additive by
    # calling is_ultrametric and is_additive with the default
    # threshold value (1e-4).
    #
    # YOUR CODE HERE
    #
    print(f"Table 1 is is_ultrametric: {is_ultrametric(dist_1)}")
    print(f"Table 2 is is_ultrametric: {is_ultrametric(dist_2)}")
    print(f"Table 1 is additive: {is_additive(dist_1)}")
    print(f"Table 2 is additive: {is_additive(dist_2)}")

    # Construct the ATPA synthase distance metric table
    atpa_table = {"1,2": 0.5, "1,3": 0.5, "1,4": 0.1, "1,5": 0.4, "1,6": 0.4,
                  "2,3": 0.3, "2,4": 0.5, "2,5": 0.5, "2,6": 0.5,
                  "3,4": 0.5, "3,5": 0.5, "3,6": 0.5,
                  "4,5": 0.4, "4,6": 0.4,
                  "5,6": 0.3}

    # Determine if the ATPA synthase distance metrics
    # are ultrametric and additive using the default
    # threshold value (1e-4).
    #
    # YOUR CODE HERE
    #
    print(f"The APTA table is ultrametric: {is_ultrametric(atpa_table)}")
    print(f"The APTA table is additive: {is_additive(atpa_table)}")


def is_ultrametric(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are ultrametric.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal
    Returns:
        (bool) True if the given distance metric is an ultrametric,
    False otherwise."""

    #
    # YOUR CODE HERE ()
    #
    # number of pairs of sequences
    dist_length = len(dist)
    # based on the number of pair of sequences to find the # of sequences
    comb_sum = 0
    count = 1
    while count < dist_length:
        comb_sum += count
        if comb_sum == dist_length:
            break
        count += 1
    # total # of sequences will be the count + 1
    num_seq = count + 1

    # loop through every possible combination of 3 sequences
    for i in range(1, num_seq-1):
        for j in range(i+1, num_seq):
            for k in range(j+1, num_seq+1):
                # (2 + 1) = 3
                str1 = str(i) + "," + str(j)
                str2 = str(i) + "," + str(k)
                str3 = str(j) + "," + str(k)
                # get the distance comparing every 2 sequences out of the 3 sequences
                d1 = dist[str1]
                d2 = dist[str2]
                d3 = dist[str3]
                # compare the 3 distances, only one pair of distances should be the same
                # if two distances are the same and one is less
                comp1 = is_almost_equal(d1, d2, threshold)
                comp2 = is_almost_equal(d1, d3, threshold)
                comp3 = is_almost_equal(d2, d3, threshold)
                count = 0
                if comp1:
                    count += 1
                if comp2:
                    count += 1
                if comp3:
                    count += 1
                # if the condition that one pair of distances should be the same holds is not true, return false
                if count != 1:
                    return False
    # every combination satisfies the condition, return true
    return True


def is_additive(dist, threshold=1e-4):
    """Check that a set of pairs of point distances are additive.

    Note: When making comparisons between distances, use `is_almost_equal` with
    the input parameterized threshold. This will be useful for subsequent
    problems where `is_ultrametric` is called. e.g. When comparing x and y, 
    also pass the threshold parameter: is_almost_equal(x, y, threshold).

    Args:
        dist (dict): exhaustive dict of pairs of points mapped to distances. 
        e.g.
            {"1,2" : 0.5, "1,3" : 0.1, "2,3" : 0.6}
        threshold (float): maximium difference in which numeric values are 
            considered equal

    Returns:
        (bool) Return True if the given distance metric is additive, 
        False otherwise."""

    #
    # YOUR CODE HERE ()
    #
    # number of pairs of sequences
    dist_length = len(dist)
    # based on the number of pair of sequences to find the # of sequences
    comb_sum = 0
    count = 1
    while count < dist_length:
        comb_sum += count
        if comb_sum == dist_length:
            break
        count += 1
    # total # of sequences will be the count + 1
    num_seq = count + 1

    # loop through every possible combination of 4 sequences
    for i in range(1, num_seq - 2):
        for j in range(i + 1, num_seq - 1):
            for k in range(j + 1, num_seq):
                for l in range(k + 1, num_seq + 1):
                    # (3 + 2 + 1) = 6
                    str1 = str(i) + "," + str(j)
                    str2 = str(i) + "," + str(k)
                    str3 = str(i) + "," + str(l)
                    str4 = str(j) + "," + str(k)
                    str5 = str(j) + "," + str(l)
                    str6 = str(k) + "," + str(l)
                    # get the distance comparing every 2 sequences out of the 4 sequences
                    d1 = dist[str1]
                    d2 = dist[str2]
                    d3 = dist[str3]
                    d4 = dist[str4]
                    d5 = dist[str5]
                    d6 = dist[str6]
                    # compare the 3 groups of distances, only one group of distance should be equal
                    # if two groups of distances are the same and one group is less
                    comp1 = is_almost_equal(d1+d6, d2+d5, threshold)
                    comp2 = is_almost_equal(d1+d6, d3+d4, threshold)
                    comp3 = is_almost_equal(d2+d5, d3+d4, threshold)
                    count = 0
                    if comp1:
                        count += 1
                    if comp2:
                        count += 1
                    if comp3:
                        count += 1
                    # if the condition that one pair of combined distances should
                    # be the same holds is not true, return false
                    if count != 1 and count != 3:
                        return False
    # every combination satisfies the condition, return true
    return True

def is_almost_equal(num_1, num_2, threshold):
    """
    Return true if the difference between the two parameters is negligible
    enough that the parameters can be considered equal.
    """
    return abs(num_1 - num_2) <= threshold


if __name__ == '__main__':
    solve_ultrametric_additive()
