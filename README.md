# computational-genomis

This repository introduces a computational perspective on the exploration and analysis of genomic and genome-scale information. It provides an integrated introduction to genome biology, algorithm design and analysis, and probabilistic and statistical modeling through 7 projects that include genome sequencing, genome sequence assembly, read mapping, local and global sequence alignment, sequence database search, gene finding, phylogenetic tree construction, and elementary gene expression analysis. Methods include dynamic programming, indexing, hidden Markov models, and elementary supervised machine learning. The projects focus on foundational algorithmic principles and development of practical experience with handling, analyzing, and visualizing genomic data using the computer language Python.

** For specific project files, check each project branch.

# Project 1: Analyze Biological Sequences

## Part 1 - Attack of the clones

We will consider an open reading frame (ORF) to be a stretch of DNA (or RNA in the case of an RNA virus like SARS-CoV-2) that starts with a start codon (ATG; or AUG in the case of RNA) and continues until it reaches the first in-frame stop codon (TAA, TAG, TGA; or UAA, UAG, UGA in the case of RNA). Specifically, we will adopt the convention that the stop codon at the end of an ORF is considered part of that ORF.

Step 1: Identification of restriction sites

Molecular cloning is a common practice in biology where bacterial cells are used to create bazillions of copies of a specific DNA fragment, often a particular gene of interest. This is accomplished by extracting the DNA fragment via polymerase chain reaction (PCR) or restriction enzyme digestion and inserting the fragment into a circular plasmid, which can then be taken up by host bacterial cells and replicated again and again. I wrote code to simulate a simplified version of this process. Specifically, cloning.py aims to locate restriction sites that flank the yeast Aim2 gene so that it can be excised from its context, identify a compatible restriction site within the multiple cloning site of the pRS304 plasmid, and then computationally insert the excised gene into the cleaved location in the plasmid.

The FASTA file aim2 plus minus 1kb.fasta contains the full genomic sequence of the yeast Aim2 gene, as well as an additional thousand base pairs upstream and downstream of the gene. This means that the part of the sequence that corresponds to the Aim2 gene itself is nucleotides 1001–1741 (inclusive).

Step 2: Cloning Aim2 into a yeast integrating plasmid

pRS304.fasta contains the sequence for a yeast integrating plasmid and pRS304 map.pdf highlights several key features of this plasmid. 

The multiple cloning site in pRS304 possesses a number of restriction sites (by design). One key step in designing a cloning experiment is determining which restriction enzyme(s) to use so that when Aim2 is excised from its flanking sequence, its sticky ends will be compatible with the sticky ends where the plasmid MCS is cut. Specifically, the overhanging single-stranded DNA sequences (often referred to as “sticky ends”) generated at the restriction sites flanking the gene need to match up with the overhangs generated at the restriction site in the plasmid.



## Part 2 - Let’s start to stop COVID-19

A critical step in stopping COVID-19 is understanding the biology of the SARS-CoV-2 virus that causes it. To do that, we can first identify all the open reading frames (ORFs) in the SARS-CoV-2 genome, demarcating them by their start and stop codons. These ORFs represent the potential proteins that the virus instructs our infected cells to make. After finding these ORFs, we can determine the amino sequences of the corresponding proteins, and continue to characterize them further. These proteins can become targets for therapies and vaccines against SARS-CoV-2 infection.

The function find_orfs in orfs.py takes as input a genome sequence and the minimum ORF length in amino acids. It returns a list of dictionaries where each dictionary entry corresponds to one ORF and contains the following information describing that ORF:
  ‘start’: the start position of the ORF (the first nucleotide of the start codon)
  ‘stop’: the stop position of the ORF (recall that we consider the stop codon part of the ORF, so the stop position will be the last nucleotide of the stop codon)
  ‘stopcodon’: specific stop codon sequence (UAG, UGA, or UAA for RNA)
  ‘nlength’: length of the ORF in nucleotides
  ‘aalength’: length of the translated peptide in amino acids
  ‘frame’: reading frame with respect to the start of the genome (0, 1, or 2)3
  ‘strand’: the strand on which the ORF is found (W or C)4. Note: For single-stranded sequences, ‘strand’ will always be set to W. However, the strand key will become more relevant in Problem 3 when we ask you to apply the ORF finder to a genomic sequence that is double-stranded.
  
Apply functions in orfs.py to the SARS-CoV-2 genome in sars cov2 wu.fasta.



## Part 3 - Plasmid Aseemly 

Step 1: The genome assembly problem

When the human genome was sequenced, researchers didn’t put it into a machine and wait for a sequence of 3 billion base pairs to come out. DNA sequencing technology (of the Sanger sequence variety, which we will discuss in more detail later) only permits the determination of around 500 to 800 base pairs at a time, and these short stretches of determined sequence are called “reads”. But we can collect these reads from bazillions of random locations in the genome. Then, to solve what is known as the genome assembly problem, algorithms are used to computationally scan for reads whose sequences partially overlap those of other reads, and then computationally merge the sequences of the overlapping reads into longer sequences of contiguous bases, nicknamed ‘contigs’ (which exist contiguously in the genome, but collectively may not include all of the genome based on which reads were collected).

I wrote plasmid.py to solve an assembly problem that is much simpler than what occurs in reality, but with one interesting twist: the DNA sequence is a circular bacterial plasmid. plasmid.fasta contains fourteen partially overlapping reads from a circular bacterial plasmid. Though the length of each read is different, their average length is approximately 625 base pairs. To simplify the problem, one can assume that each read is from the same strand of the DNA (you need not worry about reverse complementation or opposite orientations when looking for overlaps), that the end of each read overlaps the start of another read by exactly 15 base pairs (you need not worry about sliding one read against the others in all possible positions), and that there are no sequencing errors whatsoever (you need only consider perfect matches when looking for overlaps).

Function called simple_assembler reassembles a full circular plasmid sequence from any list of overlapping reads.

Step 2: Finding ORFs in a circular plasmid

Function called find_orfs_circular_double_stranded searches for ORFs within the previously reconstructed double-stranded plasmid


# Project 2: Dynamic Programming

## Part 1 - Alien invasion!

While studying an organic sample taken from an asteroid, you isolate a nucleotide sequence X from what
seems to be a mysterious alien organism. Bits of it seem to be similar to human DNA while other bits are
quite dissimilar. To better understand what kind of threat the human race might be facing, scientists around
the planet are rushing to your side to study this sequence.

One scientist, Dr. Jean Hunter, has managed to use an ORF-finder to find all the genes in X. After
translating those ORFs into peptide sequences, she has assessed, for each one, how similar it is to its closest
human counterpart, reporting a number where a more positive value indicates greater similarity and a more
negative value indicates greater dissimilarity. She has compiled all these values into one long list, where the
first element represents the human similarity for the first gene along the alien genome, the second element
represents the human similarity for the second gene along the alien genome, and so on. A very short sample
list might look like this:
                    
    [2,−3,−4, 4, 8,−2,−1, 1, 10,−5]

though the actual list for all the genes in the alien genome is quite a bit longer (owing to the fact that these
aliens may possess super-intelligence). Dr. Hunter has worked all through the night to produce this list
and is about to crash from lack of sleep, so she tasks you with the job of identifying the region (contiguous
segment) of the alien genome that possesses the highest collective similarity to the human genome. More
precisely, she defines the “collective similarity” of a region of the alien genome to be the sum of all the
similarity values of the genes contained within that region; your mission (should you choose to accept it) is
to maximize this score. For example, in the sample list above, the region of highest collective similarity is
the sublist [4, 8,−2,−1, 1, 10], whose collective similarity is 20. Note that you should return an empty region
(whose collective similarity is 0) if that ends up being optimal.

Within collective similarity.py, three different functions are implemented: brute_force, divide_and_conquer, and linear (which utilizes dynamic programming).

## Part 2 - Pebble beach

Imagine that you are given a game board arranged as a grid with n rows and 4 columns. You are also given
a set of 2n pebbles. Each pebble may be placed on at most one square of the grid, and at most one pebble
may be placed on each square of the grid; you need not use all the pebbles. Let us define a valid placement
to be a placement of some or all of the pebbles on the board such that no two pebbles lie on horizontally
or vertically adjacent squares (diagonal adjacencies are permitted). On each square of the grid is written a
positive integer. Let us define the value of a placement to be the sum of all integers written on the squares
where pebbles have been placed.

Let us consider the overall problem of finding the maximum value of a valid placement of pebbles on the
board. We say that two patterns are compatible (with one another) if they can be placed on adjacent rows
when forming a valid placement. It follows that a valid placement will consist of a succession of compatible
patterns moving down the grid. The question becomes: Among all of those valid placements, what is the
highest possible value you can get?

We can first try to break this big problem down into smaller problems by considering subproblems of size
k ∈ {0, . . . , n}, in which we consider only the first k rows of the grid. While breaking the problem down this
way is a good step, it turns out that it will not be quite enough to solve the overall problem. It will very
much help your intuition to try it for a few minutes on a small problem and see why this doesn’t quite work.
As it happens, to solve this problem we will need to consider multiple subproblems of the same size k that
differ according to their pattern type, which is the pattern that occurs in the last (kth) row. Putting this
all together, the set of subproblems we will need to solve can now be “named” based on both their size and
their pattern type.

pebble.py takes an n × 4 grid of positive integers as input and then output the maximum value of a valid placement.

The file grid.txt is a sample grid with random integers between 1 and 120. grid.txt is inputted into pebble.py for testing.


# Project 3: Genome Assembly, Short-read Mapping, 
