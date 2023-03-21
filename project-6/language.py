import random
import textwrap
import tidy_text


def run_language():
    """Call `generate_mm_text` with the provided parameters."""

    file_name = "tidy.paradise.lost.txt"
    order = 4
    M = 1000
    generated_text = generate_mm_text(file_name, order, M)
    generated_text = textwrap.wrap(generated_text, 72)

    # output the file
    filename = "paradise.lost.mm.4.txt"
    f = open(filename, "w")
    for line in generated_text:
        # format the file in lines
        f.write(line + " \n")
    f.close()
    print(generated_text)


def generate_mm_text(file_name, order, M):
    """Create a Markov model for a given text file and output artificially
    generated text from the model.

    Args:
        file_name (str): path of the text to process
        order (int): order of the Markov model
        M (int): the length of the number of characters in the returned generated
        string

    Returns:
        A string of randomly generated text using a Markov model
    """
    # Read the contents of the file
    f = open(file_name, "r")

    if f is None:
        print("Can't open " + file_name)
    else:
        contents = f.read()
        f.close()
        contents = contents.replace("\n", "")
        contents = contents.replace("\r", "")

    # Collect the counts necessary to estimate transition probabilities
    # This dictionary will store all the data needed to estimate the Markov model:
    txt_dict = collect_counts(contents, order)

    # Generate artificial text from the trained model
    seed = contents[0:order]
    text = seed

    for _ in range(M):
        next_character = generate_next_character(seed, txt_dict)
        text += next_character
        seed = seed[1:] + next_character

    text_list = textwrap.wrap(text, 72)
    text = "\n".join(text_list)

    # Return the generated text
    return text


def display_dict(txt_dict):
    """Print the text dictionary as a table of keys to counts.
    Currently accepts a dictionary specified by the return documentation in the
    `build_dict` function.

    You will need to modify this function to accept the dictionary returned by
    the `collect_counts` function.

    Arguments:
        txt_dict (dict) - Mapping keys (as strings) to counts (as ints). After
        modification for `collect_counts`, the txt_dict will map keys (as strings)
        to dictionaries of counts and followers described in the return method
        of `collect_counts`.
    """
    print("key\tcount\tfollower counts")
    for key in txt_dict.keys():
        followers_dict = txt_dict[key]['followers']
        followers_str = ""
        for follower in followers_dict:
            followers_str = followers_str + str(follower) + ":" + str(followers_dict[follower]) + " "
        print("%s\t%d\t%s " % (key, txt_dict[key]['count'], followers_str), "\t", end=" ")


def build_dict(contents, k):
    """Builds a dictionary of k-character (k-tuple) substring counts. Store the
    dictionary mapping from the k-tuple to an integer count.

    Args:
        contents (str): the string contents of to count
        k (int): number of characters in the substring

    Returns:
        a text dictionary mapping k-tuple to an integer
        Example return value with k=2:
        { 
            "ac": 1,
            "cg": 2,
            ... 
        }
    """
    #
    # YOUR CODE HERE
    #
    return_dict = {}
    # range is set this way so that the very last k-tuple in the text is not counted is the dictionary
    for i in range(len(contents)-k):
        # extract the key input
        key_input = contents[i:i+k]
        # increment by 1 if the key already exists
        if key_input in return_dict:
            return_dict[key_input] += 1
        # set value to 1 for 1st appearance if it's the first time a key appears
        else:
            return_dict[key_input] = 1
    return return_dict


def collect_counts(contents, k):
    """Build a k-tuple dictionary mapping from k-tuple to a dictionary of
    of counts and dictionary of follower counts.
    
    Args:
        contents (str): the string contents of to count
        k (int): number of characters in the substring

    Returns:
        a dictionary mapping k-tuple to a dictionary of counts and dictionary
        of follower counts. Example return value with k=2:
        {
            "ac": {
                "count": 1,
                "followers": {"g": 1, "c": 2}
            },
            ...
        }

    Note: This function will similar to `build_dict`. We separated the 
    k-character and follower counting to explain each as distinct concepts. You
    should use the k-character counting code you wrote in `build_dict` as a 
    starting point.

    While the Markov model only needs to use `collect_counts` to generate text,
    our auto-grader will test the behavior of `build_dict` so that function 
    does need to work properly.
    """
    #
    # YOUR CODE HERE
    #
    return_dict = {}
    # range is set this way so that the very last k-tuple in the text is not counted is the dictionary
    for i in range(len(contents) - k):
        # extract the key input and the next key input
        key_input = contents[i:i + k]
        key_next = contents[i + k]
        # two if conditions to account for count and followers
        if key_input in return_dict:
            return_dict[key_input]['count'] += 1
            # check if the value is already in the dictionary
            if key_next in return_dict[key_input]['followers']:
                return_dict[key_input]['followers'][key_next] += 1
            else:
                return_dict[key_input]['followers'][key_next] = 1
        # if it's a new key input, then create a new key value for the key input
        else:
            return_dict[key_input] = {}
            return_dict[key_input]['count'] = 1
            return_dict[key_input]['followers'] = {}
            return_dict[key_input]['followers'][key_next] = 1
    return return_dict


def generate_next_character(seed, txt_dict):
    """Randomly select the next character of a k-tuple using the follower
    counts to determine the probability.

    Args:
        seed (str): k-tuple to follow from
        txt_dict (dict): k-tuple count follower dictionary

    Returns:
        (str) of the next character
    """
    #
    # YOUR CODE HERE
    #

    # two lists: one for the input data structure, one for the weight of each value in the next_character_list
    next_character_list = []
    next_character_weight = []
    # extract values from the dictionaries to build the two lists
    for next_key in txt_dict[seed]['followers']:
        next_character_list.append(next_key)
        next_character_weight.append(txt_dict[seed]['followers'][next_key])

    # convert to tuple
    next_character_weight = tuple(next_character_weight)

    # use a randomized process to generate the next character based on the probability
    next_char = random.choices(next_character_list, weights=next_character_weight, k=1)
    return next_char[0]


if __name__ == "__main__":
    """Main method call, do not modify"""
    run_language()

