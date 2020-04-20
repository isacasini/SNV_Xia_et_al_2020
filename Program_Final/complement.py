# This file contains the following function(s): dna_complement

def dna_complement(sequence):
    """reads in a string of the top strand DNA sequence and returns the string of the complement (bottom strand)
    also in the 5’->3’ direction. """
    conversion = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    compsequence = "".join([conversion[c] for c in sequence][::-1])
    return compsequence