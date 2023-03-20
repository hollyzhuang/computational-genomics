from math import sqrt


def run_k_nearest_neighbor():
    """Run K Nearest neighbors against the training and test sets."""
    input_patients = read_training_data("gene_expression_training_set.txt")
    new_patients = read_test_data("gene_expression_test_set.txt")
    solve_k_output = solve_k_nearest(input_patients, new_patients, k=5)
    print(solve_k_output)
    # report the output
    for i in range(len(solve_k_output)):
        if solve_k_output[i] == "R":
            solve_k_output[i] = "Responsive to the treatment "
        if solve_k_output[i] == "N":
            solve_k_output[i] = "Non-Responsive to the treatment "
    for j in range(len(solve_k_output)):
        print(f"The prognosis for patient {j+1} is {solve_k_output[j]}")

def solve_k_nearest(input_patients, new_patients, k):
    """
    Read in the input patients as training data to use to determine and report
    the prognosis of the 10 new patients.

    Args:
        input_patients (list of dicts): dictionaries of training patient
            data. see: `read_training_data`

        new_patients (list of dicts): dictionaries of test patient data.
            see: `read_test_data`
        k (int): number of nearest neighbors to use to predict the prognosis
        for an unlabeled patient
    """
    
    #
    # YOUR CODE HERE
    #
    prognosis_list = []
    prognosis = ''
    # loop through every patient to consider each scenario
    for test in new_patients:
        # create a list to record the distance between test point and training points
        distance_list = []
        # loop through every training point
        for train in input_patients:
            # calculate the distance between each test point and training point
            d = compute_dist(test['expression'], train['expression'])
            class_name = train['class']
            # store the class name and distance for the 1000 points
            distance = {"class": class_name, "distance": d}
            distance_list.append(distance)
        # sort the list from small to large based on the distance value
        distance_list = sorted(distance_list, key=lambda x: x['distance'])
        # get the highest k values
        distance_list = distance_list[:k]
        # check which one is the majority, response or non-responsive
        count_r = 0
        count_n = 0
        for i in distance_list:
            if i['class'] == 'R':
                count_r += 1
            if i['class'] == 'N':
                count_n += 1
        # responsive wins
        if count_r > count_n:
            prognosis = 'R'
        # non-responsive wins
        if count_n > count_r:
            prognosis = 'N'
        prognosis_list.append(prognosis)
    return prognosis_list


def read_training_data(file_name):
    """Read the training gene expression data from a text file. Note: the
    patients in the training data are classified as "R" (responsive to
    treatment) or "N" (non-responsive to treatment).  For example,
    input_patients[0]["class"] = the class of the first patient (R or N)
    input_patients[0]["expression"][0] = the expression of the first
    gene for the first patient.

    Returns:
        (list of dicts): list of patients as a class and expression data. The
        dictionary of each patient will be in the form of:
            'class' -> string with values strictly 'N' or 'R' for
            non-responsive or responsive to the treatment
            'expression' -> list of floats of gene expression values

        and look something like:
            {'class': 'N', 'expression': [9.049, 8.313, ..., 6.428700888]}
    """
    return read_data(file_name, test_data=False)


def read_test_data(file_name):
    """Read the test gene expression data from a text file. Note: the
    patients in the test data are not classified.

   Returns:
    (list of dicts): list of patients as a class and expression data. The
    dictionary of each patient will be in the form of:
        'class' -> string with only 'unknown' as its value
        'expression' -> list of floats of gene expression values

    and look something like:
        {'class': 'unknown', 'expression': [9.049, 8.313, ..., 6.428700888]}
    """
    return read_data(file_name, test_data=True)


def read_data(file_name, test_data=False):
    with open(file_name, "r") as f:
        lines = f.readlines()

    patients = []

    for line in lines:
        line = line.strip()
        data = line.split()  # check that you are splitting on "\t"

        if test_data:
            class_name = "unknown"
        else:
            class_name = data.pop(0)

        float_data = [float(datum) for datum in data]
        patient = {"class": class_name, "expression": float_data}
        patients.append(patient)

    return patients


def compute_dist(tuple_1, tuple_2):
    """Return the Euclidean distance between two points in any number of
    dimensions."""
    
    if len(tuple_1) != len(tuple_2):
        raise ValueError("Cannot compute Euclidean distance between tuples of different sizes!")
    
    dist = 0
    for i in range(len(tuple_1)):
        dist += (tuple_1[i] - tuple_2[i]) * (tuple_1[i] - tuple_2[i])
      
    return sqrt(dist)


if __name__ == "__main__":
    run_k_nearest_neighbor()
