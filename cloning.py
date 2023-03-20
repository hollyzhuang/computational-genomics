import re
from typing import Union, Any

from compsci260lib import get_fasta_dict


def run_cloning():
    """Function demonstrating the procedure of extracting a gene and inserting
    into a plasmid using restriction enzymes.
    """
    aim2_fasta_path = "aim2_plus_minus_1kb.fasta"
    pRS304_fasta_path = "pRS304.fasta"

    # Read in the Aim2 genomic sequence from the fasta file along with its
    # upstream and downstream regions.
    #
    # Your code goes here
    #
    aim2_dict = get_fasta_dict(aim2_fasta_path)
    aim2_genomic = aim2_dict["aim2"].upper()

    # Store the beginning and end locations of the Aim2 gene as Python indices.
    aim2_beg = 1001 - 1
    aim2_end = 1741 - 1

    # Report start, end, length in nucleotides, and length in amino acids.
    #
    # Your code goes here
    #
    # nucleotide length
    aim2_length_nucleotides = aim2_end - aim2_beg + 1
    # amino acid length
    aim2_length_aa = int((aim2_length_nucleotides - 3) / 3)

    print(aim2_beg, aim2_end, aim2_length_nucleotides, aim2_length_aa)
    print("The starting and ending locations for the gene are respectively " + str(aim2_beg + 1) + " and " + str(aim2_end + 1) + ".")
    print("The length of gene in nucleotides is " + str(aim2_length_nucleotides) + ".")
    print("The length of gene in amino acids is " + str(aim2_length_aa) + ".")

    # Define regular expression terms for each restriction enzyme
    r_enzymes = get_restriction_enzymes_regex()

    # Store coordinates of restriction sites found upstream, downstream, and
    # within the aim2 gene
    r_enzyme_sites = find_restriction_sites(aim2_genomic, window_beg=aim2_beg,
                                            window_end=aim2_end, r_enzymes=r_enzymes)

    # Report the found restriction enzyme sites
    #
    # Your code goes here
    #
    print(r_enzyme_sites)

    # Read in the pRS304 plasmid sequence and find restriction sites
    # in the entire plasmid
    mcs_start = 1891 - 1
    mcs_end = 1993 - 1
    prs304_genomic = (get_fasta_dict(pRS304_fasta_path))["pRS304"].upper()
    all_prs304_sites = find_restriction_sites(prs304_genomic, window_beg=mcs_start,
                                              window_end=mcs_end, r_enzymes=r_enzymes)
    print(all_prs304_sites)
    # Report relevant summary information
    summary = report_pRS304_MCS_sites(all_prs304_sites)
    print(summary)


    # Extract aim2 gene and insert into the plasmid, report its length
    #
    # Your code goes here
    #
    aim2_start_cut = 458 + 1
    aim2_end_cut = 2611 - 1
    aim2_bp_start = aim2_start_cut + 1
    aim2_bp_end = aim2_end_cut + 1
    aim2_cut = aim2_genomic[aim2_start_cut:aim2_end_cut+1]
    print("The nucleotide positions (with respect to the original sequence) for the beginning and end of the excised fragment are respectively " + str(aim2_bp_start) + " and " + str(aim2_bp_end))
    print(aim2_cut)

    prs304_break_start = 1922
    prs304_break_end = 1927
    prs304_genomic_1 = prs304_genomic[0:prs304_break_start + 1]
    prs304_genomic_2 = prs304_genomic[prs304_break_end:len(prs304_genomic)]
    inserted = prs304_genomic_1 + aim2_cut + prs304_genomic_2
    length = len(inserted)
    print(length)


def get_restriction_enzymes_regex():
    """Returns a dictionary of restriction enzyme regular expressions for
    searching in genomic sequences.

    This function should be used for find_restriction_sites.
    """

    return {
        "BamHI": "GGATCC",
        "BstYI": "rGATCy",
        "SpeI": "ACTAGT",
        "SphI": "GCATGC",
        "StyI": "CCwwGG",
        "SalI": "GTCGAC"
    }


def find_restriction_sites(genomic_seq, window_beg, window_end, r_enzymes):
    """Finds the restriction enzyme sites in a genomic sequence. Stored
    as a dictionary of lists of dictionaries. Each restriction enzyme is
    a key in the outer dictionary corresponding to a list of all its sites.
    Each site has its information represented in a dictionary. Each site
    entry has a key specifying if the site is downstream, upstream, or
    within the parameter-specified window.

    Args:
        genomic_seq (str): genomic sequence to search for sites in
        window_beg (int)
        window_end (int)
        r_enzymes (dict): dictionary where key:value is enzyme:regex (same
        format returned by the get_restriction_enzymes_regex function)

    Each found site will have defined keys:
        sequence (str): the nucleotide sequence matched in the genome
        start (int): the start nucleotide position in the genome
        end (int): the ending position
        location (str): the position of the site relative to the specified window.
            (must be "upstream", "downstream", or "within")

    A valid returned dictionary may look like this:
    {
        "BamHI" : [
            {
                "start": 10,
                "end": 15,
                "sequence": "GGATCC",
                "location": "upstream"
            },
            {
                "start": 100,
                "end": 105,
                "sequence": "GGATCC",
                "location": "downstream"
            }
        ],
        "BstYI" : [
            {
                "start": 30,
                "end": 35,
                "sequence": "AGATCC",
                "location": "within"
            }
        ]
    }
    """
    if type(genomic_seq) is not str:
        raise TypeError(f"the argument 'genomic_seq' must be a string. "
                        "(received {type(genomic_seq)})")

    if window_beg >= window_end:
        raise ValueError("window_beg must be smaller than window_end")

    # Create dictionary to store coordinates
    r_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}

    #
    # Your code goes here
    #
    # loop through every restriction enzyme
    for key in r_enzymes:
        # loop through the genomic sequence
        for i in range(len(genomic_seq) - len(r_enzymes[key])):
            # special case for BstYI
            if "r" in r_enzymes[key]:
                # location 1 and location 6 of the restriction enzyme can have diff possibilities
                ry = genomic_seq[i] + genomic_seq[i+5]
                # list all the possible combinations
                ry_possible = ["AC", "AT", "GC", "GT"]
                # look for matching sequence, taking ry_possible into account
                if genomic_seq[(i + 1):(i+len(r_enzymes[key]) - 1)] == (r_enzymes[key])[1:5] and ry in ry_possible:
                    # determine if the location is upstream, within, or downstream based on
                    # comparison with the start and end locations
                    if i < window_beg:
                        location = "upstream"
                    elif i <= window_end:
                        location = "within"
                    else:
                        location = "downstream"
                    # create new dictionary with stored values and append it to the return enzyme dictionary
                    new_dict = {"start": i + 1, "end": i + 6, "sequence": genomic_seq[i: i+6], "location": location}
                    r_enzyme_sites[key].append(new_dict)
            # special case for StyI
            if "w" in r_enzymes[key]:
                # location 3 and 4 of the restriction enzyme can have diff possibilities
                ww = genomic_seq[i+2:i+4]
                # list all the possible combinations
                ww_possible = ["AA", "TT", "AT", "TA"]
                # look for matching sequence, taking ww_possible into account
                if genomic_seq[i:i + 2] == (r_enzymes[key])[0:2] and genomic_seq[i + 4:i+len(r_enzymes[key])] == (r_enzymes[key])[4:6] and ww in ww_possible:
                    # determine if the location is upstream, within, or downstream based on
                    # comparison with the start and end locations
                    if i < window_beg:
                        location = "upstream"
                    elif i <= window_end:
                        location = "within"
                    else:
                        location = "downstream"
                    # create new dictionary with stored values and append it to the return enzyme dictionary
                    new_dict = {"start": i + 1, "end": i + 6, "sequence": genomic_seq[i: i+6], "location": location}
                    r_enzyme_sites[key].append(new_dict)
            # the rest of the case - look for matching restriction enzyme sequence
            if genomic_seq[i:i+len(r_enzymes[key])] == r_enzymes[key]:
                # determine if the location is upstream, within, or downstream based on
                # comparison with the start and end locations
                if i < window_beg:
                    location = "upstream"
                elif i <= window_end:
                    location = "within"
                else:
                    location = "downstream"
                # create new dictionary with stored values and append it to the return enzyme dictionary
                new_dict = {"start": i + 1, "end": i + 6, "sequence": genomic_seq[i: i+6], "location": location}
                r_enzyme_sites[key].append(new_dict)

    return r_enzyme_sites


def report_pRS304_MCS_sites(p_enzyme_sites):
    """For each restriction enzyme, report how often that enzyme cuts the plasmid
    outside the MCS (upstream or downstream), how often it cuts the plasmid inside
    (within) the MCS, and relevant details about any sites located inside the MCS.
    """
    #
    # Your code goes here
    #
    # create a new dictionary to store the information for each restriction enzyme
    mcs_enzyme_sites = {"BamHI": [], "BstYI": [], "SpeI": [],
                      "SphI": [], "StyI": [], "SalI": []}
    # loop through the dictionary with the cut sites for each restriction enzyme
    for p in p_enzyme_sites:
        count_in = 0
        count_out = 0
        detail_list = []
        # loop through every cut site of the specific restriction enzyme
        for i in range(len(p_enzyme_sites[p])):
            cut_site = (p_enzyme_sites[p])[i]
            location = cut_site["location"]
            # check to see if the location is within or upstream/downstream,
            # then count their respective frequency
            if location == "within":
                count_in += 1
                # append to a list if cut inside mcs
                detail_list.append(cut_site)
            else:
                count_out += 1
        # create a new dictionary to record cut sites that are located outside mcs, inside mcs,
        # and the detailed information of sites cut within mcs
        new_dict = {"outside_mcs": count_out, "inside_mcs": count_in, "details": detail_list}
        mcs_enzyme_sites[p].append(new_dict)

    return mcs_enzyme_sites


if __name__ == "__main__":
    """Run run_cloning(). Do not modify this code"""
    run_cloning()
