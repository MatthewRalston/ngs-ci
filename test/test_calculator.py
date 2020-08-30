import os
import sys
import unittest

sys.path.append("..")

from functools import reduce
import pysam
from ngsci import calculator, parser

class Test_dissimilarity(unittest.TestCase):
    def setUp(self):
        self.chrom = "NC_001988.2"
        self.samfile = pysam.AlignmentFile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/test.bam"), "rb")
        iter = self.samfile.fetch(self.chrom, 0, 10)
        self.reads = [x for x in iter]
        self.read1 = self.reads[0]
        self.read2 = self.reads[2]

    def tearDown(self):
        self.samfile.close()

    def test_throws_TypeError(self):
        """
        ngsci.calculator.dissimilarity throws a TypeError when one or both reads aren't pysam.AlignedSegment objects
        """
        with self.assertRaises(TypeError):
            calculator.dissimilarity(None, self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity("hello", self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(1, self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(1.0, self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity([], self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity({}, self.read2)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, None)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, "hello")
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, 1)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, 1.0)
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, [])
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, {})

        # Also throws TypeError if both reads are the same
        with self.assertRaises(TypeError):
            calculator.dissimilarity(self.read1, self.read1)

    def test_returns_an_integer(self):
        """
        The dissimilarity function returns an integer value
        """
        self.assertEqual(type(calculator.dissimilarity(self.read1, self.read2)), int)

    def test_calculates_dissimilarity(self):
        """
        Calculates the pairwise dissimilarity (unique bases) between two reads
        """
        self.assertEqual(calculator.dissimilarity(self.read1, self.read2), 3)


def tri(x, n=0):
    if x == 0:
        return n
    else:
        return tri(x-1, n+x)
        

class Test_max_summed_dissimilarity(unittest.TestCase):
    def setUp(self):
        self.read_length = 76

    def test_throws_TypeError(self):
        """
        ngsci.calculator.max_summed_dissimilarity throws a TypeError when the read_length is not an int
        """
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity(None)
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity("hello")
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity(True)
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity(1.0)
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity([])
        with self.assertRaises(TypeError):
            calculator.max_summed_dissimilarity({})


    def test_is_triangular_sum(self):
        """
        This function returns a sum of triangular numbers that has an equivalent algebraic expression. This test checks that the triangular sum is equivalent to the algebraic formulation used by the function.
        """
        triangular_sum = reduce(lambda x,y: x+y, map(lambda z: tri(self.read_length - z) + tri(z - 1), range(1, self.read_length+1)))
        algebraic_max_summed_dissimilarity = calculator.max_summed_dissimilarity(self.read_length)
        # Assert that the algebraic simplification used in the calcultor is equal to a formulation based on triangular numbers
        self.assertEqual(triangular_sum, algebraic_max_summed_dissimilarity)
        #  Empirically, the max_summed_dissimilarity of 76 is 146300
        self.assertEqual(algebraic_max_summed_dissimilarity, 146300)


class Test_complexity_index(unittest.TestCase):
    def setUp(self):
        self.chrom = "NC_001988.2"
        self.samfile = pysam.AlignmentFile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/saturated.bam"), "rb")
        iter = self.samfile.fetch(self.chrom, 75, 76)
        allreads = [x for x in iter]
        self.bamreader = parser.BamReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/test.bam"))
        self.positive, self.negative = self.bamreader._partition_reads(allreads, False)

    def tearDown(self):
        self.samfile.close()
        self.bamreader.samfile.close()
        
    def test_reads_is_list(self):
        """
        ngsci.calculator.complexity_index throws a TypeError when reads is not a list
        """
        with self.assertRaises(TypeError):
            calculator.complexity_index(None, 75, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index("hello", 75, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index(1, 75, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index(1.0, 75, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index({}, 75, 10500)
    
    def test_read_length_is_int(self):
        """
        ngsci.calculator.complexity_index throws a TypeError when read_length is not an int
        """
        with self.assertRaises(TypeError):
            calculator.complexity_index([], None, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], "hello", 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 1.0, 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], [], 10500)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], {}, 10500)

    def test_denominator_is_int(self):
        """
        ngsci.calculator.complexity_index throws a TypeError when the denominator is not an int
        """
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 76, None)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 76, "hello")
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 76, 1.0)
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 76, [])
        with self.assertRaises(TypeError):
            calculator.complexity_index([], 76, {})

    def test_saturated_complexity_index_is_100(self):
        """
        complexity_index should calculate the complexity of this read group as 100%
        """
        self.assertEqual(calculator.complexity_index(self.positive, 76, self.bamreader.denominator), 100)
            
if __name__ == '__main__':
    unittest.main()
