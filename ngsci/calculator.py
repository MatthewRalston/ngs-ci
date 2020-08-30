import os
import itertools
import functools

import logging
logger = logging.getLogger(__file__)

from pysam import AlignedSegment


def dissimilarity(read1: AlignedSegment, read2: AlignedSegment):
    """
    Calculates the dissimilarity between two reads

    :param read1: The first read to be compared
    :type read1: pysam.AlignedSegment
    :param read2: The second read to be compared
    :type read2: pysam.AlignedSegment
    :returns: Length of non-overlapping/unique bases
    :rtype: int
    """
    if not isinstance(read1, AlignedSegment):
        raise TypeError("ngsci.calculator.dissimilarity expects a pysam.AlignedSegment as its first positional argument")
    elif not isinstance(read2, AlignedSegment):
        raise TypeError("ngsci.calculator.dissimilarity expects a pysam.AlignedSegment as its second positional argument")
    elif read1 == read2:
        raise TypeError("ngsci.calculator.dissimilarity expects two difference pysam.AlignedSegment reads as its arguments.")
    elif read1.is_unmapped or read2.is_unmapped:
        raise TypeError("Cannot calculate the dissimilarity if one or both reads are unmapped.")
    if read1.reference_start == read2.reference_start:
        return 0
    elif read1.reference_start > read2.reference_start: # Read2 Read1
        if read1.reference_end < read2.reference_end: # Read1 inside read2!! (Bad, unequal read lengths)
            # logger.error("         In the comparison of two reads:")
            # logger.error("Read1 start: {0}   stop: {1}".format(read1.reference_start, read1.reference_end))
            # logger.error("Read2 start: {0}   stop: {1}".format(read2.reference_start, read2.reference_end))
            # logger.error("Read 1 is inside read 2. This makes proper calculation of the dissimilarity impossible")
            # raise IndexError("ngsci.calculator.dissimilarity expects reads to have uniform read lengths.")
            return (read1.reference_start - read2.reference_start) + (read2.reference_end - read1.reference_end)
        else: # Normal overlap
            return read1.reference_start - read2.reference_start
    else: # Read1 Read2 <- appearance on the genome
        if read1.reference_end > read2.reference_end:
            # logger.error("         In the comparison of two reads:")
            # logger.error("Read1 start: {0}   stop: {1}".format(read1.reference_start, read1.reference_end))
            # logger.error("Read2 start: {0}   stop: {1}".format(read2.reference_start, read2.reference_end))
            # logger.error("Read 2 is inside read 1. This makes proper calculation of the dissimilarity impossible")
            # raise IndexError("ngsci.calculator.dissimilarity expects reads to have uniform read lengths.")
            return (read2.reference_start - read1.reference_start) + (read1.reference_end - read2.reference_end)
        else:
            return read2.reference_start - read1.reference_start

def summed_dissimilarity(reads: list):
    """
    Calculates the summed dissimilarity across a group of reads

    :param reads: A list of pysam.AlignedSegment reads to calculate dissimilarities
    :type reads: list
    :returns:
    :rtype: int
    """
    if type(reads) is not list:
        raise TypeError("ngsci.calculator.summed_dissimilarity expects a list of reads as it first positional argument")
    elif not all(isinstance(x, AlignedSegment) for x in reads):
        raise TypeError("ngsci.calculator.summed_dissimilarity expects a list of pysam.AlignedSegment as its first positional argument")
    numreads = len(reads)
    
    sum = 0
    for x, y in itertools.permutations(reads, 2):
        logger.debug("Read1: {0}\nRead2: {1}".format(x, y))
        sum += dissimilarity(x, y)
    return sum


def max_summed_dissimilarity(read_length: int):
    """Calculates the max_summed_dissimilarity for a given read_length. The algebraic formulation is a little complicated, and the triangular number equivalent is equally complicated. I'd uggest reading the whitepaper for the full formulation.

    :param read_length: The length of reads of the sequencing library.
    :type read_length: int
    :returns: maximum possible summed dissimilarity
    :rtype: int
    """
    if type(read_length) is not int:
        raise TypeError("ngsci.calculator.max_summed_dissimilarity expects an int as its first positional argument")
    return int(functools.reduce(lambda x,y: x+y, map(lambda r: (read_length**2)/2 - (read_length*r) + (read_length/2) + (r**2) - r, range(1, read_length+1))))
    

def denominator_calc(read_length: int):
    """Calculates the denominator of the complexity index as the product of the maximum_summed_dissimilarity of that read_length and the read_length itself.

    :param read_length: The length of reads of the sequencing library.
    :type read_length: int
    :returns: the denominator including normalization factors for the complexity index
    :rtype: int
    """
    return read_length*max_summed_dissimilarity(read_length)

def complexity_index(reads: list, read_length: int, denominator: int):
    """Calculates the fully-normalized complexity index for a list of deduplicated/unique reads

    :param reads: A list of deduplicated pysam.AlignedSegment objects
    :type reads: list
    :param read_length: The read_length sampled from the BAM file.
    :type read_length: int
    :param denominator: The denominator estimated from the read_length
    :type denominator: int
    :returns: the complexity index of the group of reads.
    :rtype: int
    """
    if type(reads) is not list:
        raise TypeError("ngsci.calculator.complexity_index expects a list as its first positional argument")
    elif type(read_length) is not int:
        raise TypeError("ngsci.calculator.complexity_index expects an int as its second positional argument")
    elif type(denominator) is not int:
        raise TypeError("ngsci.calculator.complexity_index expects an int as its third positional argument")
    return int(round(100*len(reads)*summed_dissimilarity(reads)/denominator))
