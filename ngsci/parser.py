import os
import sys

import logging
logger = logging.getLogger(__file__)
#logger.setLevel(logging.DEBUG)

import numpy as np
import pysam
from ngsci import calculator

class BamReader:
    """
    :ivar bamfile: The SAM/BAM/CRAM file to process
    :ivar strand_specific: The F/FR/RF strand_specific chemistry
    """

    def __init__(self, bamfile: str):
        if type(bamfile) is not str:
            raise TypeError("ngsci.reads.BamReader expects a str as its first positional argument")
        elif not os.path.exists(bamfile):
            raise IOError("ngsci.reads.BamReader expects an existing SAM/BAM/CRAM file as its first positional argument")

        # Calculate read length using samtools to subsample
        self.read_length = None
        self.denominator = None
        for x in [".00001", ".0001", ".001", ".01", ".1"]:
            reads = pysam.view("-s", "1234{0}".format(x), bamfile)
            if reads != '':
                self.read_length = max(set(map(lambda r: len(r.split("\t")[9]), reads.rstrip("\n").split("\n"))))
                break
        if self.read_length is None:
            raise RuntimeError("ngsci.reads.BamReader could not determine read length when sub-sampling as much as 10% of the file")
        else:
            logger.warning("Read length determined to be: {0}".format(self.read_length))
            self.denominator = calculator.denominator_calc(self.read_length)

        suffix = os.path.splitext(bamfile)[-1]
        if suffix == ".sam":
            self.samfile = pysam.AlignmentFile(bamfile, "r")
        elif suffix == ".bam":
            self.samfile = pysam.AlignmentFile(bamfile, "rb")
        elif suffix == ".cram":
            self.samfile = pysam.AlignmentFile(bamfile, "rc")
        else:
            raise IOError("ngsci.reads.BamReader expects a .sam, .bam, or .cram file as its first positional argument")
        try:
            if not self.samfile.check_index():
                raise IOError("ngsci.reads.BamReader requires an indexed bam/cram file as its first positional argument")
        except ValueError as e:
            logger.error(e)
            raise IOError("ngsci.reads.BamReader requires an indexed bam/cram file as its first positional argument")
        self.chroms = set(map(lambda x: x.contig, self.samfile.get_index_statistics()))

            
    def fetch(self, chrom:str, start:int, step:int, strand_specific=False):
        """
        Calculates the complexity index across a region specified by a start site and 'step' bases in the 5'->3' direction on the forward strand. The complexity index can be calculated in a strand-specific manner.

        :param chrom: The chromosome to calculate complexities from.
        :type chrom: str
        :param start: The start site in the chromosome to fetch from.
        :type start: int
        :param step: The step size or number of complexities to calculate downstream from the start site.
        :type step: int
        :param strand_specific: Whether or not to calculate complexities in a strand-specific manner.
        :type strand_specific: bool
        :returns: a 2-membered tuple of complexity indices for the positive and negative strands, respectively
        :rtype: tuple
        """
        if type(chrom) is not str:
            raise TypeError("ngsci.reads.BamReader.fetch expects a str as its first positional argument")
        elif chrom not in self.chroms:
            logger.error("Chromosome '{0}' not in the list of chromosomes:\n{1}".format(chrom, self.chroms))
            raise TypeError("ngsci.read.BamReader.fetch expects a valid contig identifier as its first positional argument")
        elif type(start) is not int or start < 0:
            raise TypeError("ngsci.reads.BamReader.fetch expects a positive int as its second positional argument")
        elif type(step) is not int or step < 0:
            raise TypeError("ngsci.reads.BamReader.fetch expects a positive int as its second positional argument")
        elif type(strand_specific) is not bool:
            raise TypeError("ngsci.reads.BamReader.fetch expects the keyword argument 'strand_specific' to be a bool")
        iter = self.samfile.pileup(chrom, start, start+step, truncated=False, stepper="nofilter", ignore_overlaps=False)
        positive = np.zeros(step).astype('uint16')
        negative = np.zeros(step).astype('uint16')
        i = 0
        for col in iter:
            logger.debug("-- Fetching reads from position '{0}'".format(col.reference_pos))
            if i < step:
                pos, neg = self._get_complexity(list(map(lambda x: x.alignment, col.pileups)), strand_specific)
                positive[i] = pos
                negative[i] = neg
            i += 1
        return positive, negative


    
            
    def _get_complexity(self, reads:list, strand_specific:bool):
        """
        Calculates complexity of all reads aligned to a specific base

        :param reads: A list of reads aligned to a specific base
        :type reads: list
        :params strand_specific: whether or not to calculate strand-specific complexity
        :type strand_specific: bool
        :returns: The integers of complexity for the +/- strands, respectively
        :rtype: tuple
        """
        positive, negative = self._partition_reads(reads, strand_specific)
        return calculator.complexity_index(positive, self.read_length, self.denominator), calculator.complexity_index(negative, self.read_length, self.denominator)
        
    def _partition_reads(self, reads: list, strand_specific: bool):
        """
        Partitions and unique-ifys reads to positive and negative strands (optional)


        :param reads: A list of reads aligned to a specific base
        :type reads: list
        :param strand_specific: whether or not to partition reads according to +/- strand
        :type strand_specific: bool
        :returns: The tuple of lists of reads aligned to the +/- strand, respectively
        :rtype: tuple
        """
        if type(reads) is not list:
            raise TypeError("ngsci.reads._partition_reads expects a list of pysam.AlignedSegment objects as its first positional argument")
        positive = []
        negative = []
        pos = set()
        neg = set()
        for r in reads:
            if r.is_unmapped:
                continue
            elif strand_specific:
                # Position hasn't been seen. Add
                if r.reference_start not in pos and r.reference_start not in neg:
                    if r.is_read1:
                        if r.is_reverse:
                            negative.append(r)
                            neg.add(r.reference_start)
                        else:
                            positive.append(r)
                            pos.add(r.reference_start)
                    else: # Is read2
                        if r.is_reverse:
                            positive.append(r)
                            pos.add(r.reference_start)
                        else:
                            negative.append(r)
                            neg.add(r.reference_start)
                # Read has been seen on pos or neg before, check strand
                elif r.reference_start in pos:
                    if r.is_read1:
                        # Read is correctly on negative strand and hasn't been seen
                        if r.is_reverse: 
                            negative.append(r)
                            neg.add(r.reference_start)
                    else: # Read2 refers to a + strand fragment in FR chemistry if it's on the negative strand
                        if not r.is_reverse:
                            negative.append(r)
                            neg.add(r.reference_start)
                elif r.reference_start in neg:
                    if r.is_read1:
                        # Read is correctly on negative strand and hasn't been seen
                        if not r.is_reverse: 
                            positive.append(r)
                            pos.add(r.reference_start)
                    else: # Read2 refers to a + strand fragment in FR chemistry if it's on the negative strand
                        if r.is_reverse:
                            positive.append(r)
                            pos.add(r.reference_start)
                else:
                    continue
            else:
                if r.reference_start not in pos:
                    positive.append(r)
                    pos.add(r.reference_start)
                else: # Reference position already seen, move on
                    continue
        logger.debug("    From {0} reads, {1} unique on positive strand, {2} unique on negative.".format(len(reads), len(positive), len(negative)))
        return positive, negative

