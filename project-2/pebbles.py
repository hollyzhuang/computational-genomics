def solve_pebbles(grid_file):
    """Code for the "Pebble Beach" problem. This problem involves implementing
    an O(n) dynamic programming algorithm for computing the maximum value of
    the placement of pebbles under the constraint that no pebbles can be
    vertically or horizontally adjacent.

    Args: grid_file (str): a string with the name of the file that contains
          the grid of values. Each line of that file should contain a row
          of four integers, separated by tabs

    Returns: the maximal score for the optimal pebble placements
    """

    #
    # Your code goes here
    #
    # read the file
    f = open(grid_file, "r")
    # create an empty list for the matrix
    input_grid = []
    n = 0
    # process the file into a matrix
    for line in f:
        n += 1
        line_list = line.strip().split()
        input_grid.append(line_list)
    # the numbers are string values, so convert them to integer
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            input_grid[i][j] = int(input_grid[i][j])
    # 8 different arrangements for each row
    # Pattern 0: No peddle is placed in any of the column
    # Pattern 1-4: One pebble placed in 1, 2, 3 4 column respectively and all other columns are empty
    # Pattern 5: One pebble placed in column 1 and another placed in column 3
    # Pattern 6: One pebble placed in column 1 and another placed in column 4
    # Pattern 7: One pebble placed in column 2 and another placed in column 4
    # create an index list to specify different patterns
    pattern_index = {0: [], 1: [0], 2: [1], 3: [2], 4: [3], 5: [0, 2], 6: [0, 3], 7: [1, 3]}
    # create a list that indicates what patterns each pattern is compatible with
    compatible_list = {0: [1, 2, 3, 4, 5, 6, 7], 1: [0, 2, 3, 4, 7], 2: [0, 1, 3, 4, 5, 6],
                       3: [0, 1, 2, 4, 6, 7], 4: [0, 1, 2, 3, 5], 5: [0, 2, 4, 7],
                       6: [0, 2, 3], 7: [0, 1, 3, 5]}

    # create a matrix list to store the n by 8 max value valid placement
    matrix_list = [[]]

    # base case values for every combination of the 8 placements
    # append these 8 values to the first row of the matrix
    for i in range(8):
        value = 0
        for j in pattern_index[i]:
            value += input_grid[1][j]
        matrix_list[0].append(value)

    # traverse through 1 to n because at index 0 it's already calculated
    for i in range(1, n):
        # create a new list every time a row is iterated
        matrix_list.append([])
        # traverse through 8 different possible arrangements
        for j in range(8):
            covered_pebbles = 0
            # get the covered pebbles values at each j
            for m in pattern_index[j]:
                covered_pebbles += input_grid[i][m]
            # get the max former total sum for Pattern j placement at index i-1
            max_former_sub = 0
            for compatible in compatible_list[j]:
                max_former_sub = max(matrix_list[i - 1][compatible], max_former_sub)
            # append the value to A(i, j)
            matrix_list[i].append(max_former_sub + covered_pebbles)

    # traverse through the last row of the matrix and get the maximum value for a valid placement
    max_sum = 0
    for k in matrix_list[n-1]:
        max_sum = max(k, max_sum)

    return max_sum


def run_pebbles():
    """Run solve pebbles, you may try an create different grid files to debug
    that match the formatting of grid.txt (tab separated values)

    Note: there is no reporting in this problem, only `solve_pebbles` will
    be evaluated for correctness.
    """
    max_score = solve_pebbles('grid.txt')
    print("The max score for this grid is %d" % max_score)


if __name__ == '__main__':
    """Run run_pebbles(). Do not modify this code"""
    run_pebbles()
