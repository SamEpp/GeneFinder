# -*- coding: utf-8 -*-
"""
This is a gene finder. It will identify stings of protiens by finding
the start and end condons as well as the frame of refrence.
@author: Samantha Eppinger

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('')
    'Not valid nucleotide'
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'Not valid nucleotide'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    rev_dna = ''
    for index in range(len(dna)-1,-1,-1):
        comp = get_complement(dna[index])
        rev_dna += comp
    return rev_dna

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGAAAATAA")
    'ATGAGA'
    """
    for index in range(0,len(dna),3):
        dna_section = dna[index: index+3]
        if  dna_section == 'TAG' or dna_section == 'TGA' or dna_section == 'TAA':
            return dna[:index]
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    index = 0
    dna_start_stop_list = []
    while index < len(dna)-2:
        dna_start = dna[index: index+3]
        if dna_start == 'ATG':
            orf = rest_of_ORF(dna[index:])
            dna_start_stop_list.append(orf)
            index = index + len(orf)
        else:
            index = index + 3
    return dna_start_stop_list

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCATGAATGTAGAAGATG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    index = 0
    total_orf = []
    while index < 3:
        one_orf_list = find_all_ORFs_oneframe(dna[index:])
        total_orf = total_orf+(one_orf_list)
        index = index + 1
    return total_orf

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    both_strands_list = []
    both_strands_list = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return both_strands_list

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest_strand = ""
    for orf in find_all_ORFs_both_strands(dna):
        if len(orf) > len(longest_strand):
            longest_strand = orf
    return longest_strand
    # This will only keep track of the longest strand
    #(not keep track of the shorter lines)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest_strand = 0
    i = 0
    while i < num_trials:
        if len(longest_ORF(dna)) > longest_strand:
            longest_strand = len(longest_ORF(dna))
        dna = shuffle_string(dna)
        i = i+1
    return longest_strand

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """

    protien_string= ""
    for index in range(0, len(dna), 3):
        if len(dna[index:index+3]) == 3:
            amino_acid = aa_table[dna[index:index+3]]
            protien_string += amino_acid
    return protien_string

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    y = ""
    protien_long_list = []
    # took our the following line due to the impossibly high number it would return.
    # m = longest_ORF_noncoding(dna, 1500)
    m = 200
    x_list = find_all_ORFs_both_strands(dna)
    for x in x_list:
        if len(x) > m:
            protien_long_list.append(coding_strand_to_AA(x))
    return protien_long_list



if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    # the print command will make running the program easier
    # Tested the program in "blast" and got a 100% match for salmonella
    print(gene_finder(dna))
    """
    doctest.run_docstring_examples(coding_strand_to_AA, globals())
    """
