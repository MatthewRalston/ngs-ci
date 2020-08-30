import os
import sys
import unittest

sys.path.append("..")

from functools import reduce
import pysam
from ngsci import parser, calculator

class Test_BamReader(unittest.TestCase):
    def setUp(self):
        self.chrom = "NC_001988.2"
        self.bamreader = parser.BamReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/test.bam"))
        self.samfile = pysam.AlignmentFile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/saturated.bam"))
        

    def tearDown(self):
        self.bamreader.samfile.close()
        self.samfile.close()

    def test_BamReader_expects_a_str(self):
        """
        ngsci.reads.BamReader throws a TypeError when the input file isn't a str
        """
        with self.assertRaises(TypeError):
            parser.BamReader(None)
        with self.assertRaises(TypeError):
            parser.BamReader(1)
        with self.assertRaises(TypeError):
            parser.BamReader(1.0)
        with self.assertRaises(TypeError):
            parser.BamReader([])
        with self.assertRaises(TypeError):
            parser.BamReader({})
            
    def test_BamReader_expects_existing_file(self):
        """
        ngsci.reads.BamReader throws an IOError when the str isn't an existing filepath
        """
        with self.assertRaises(IOError):
            parser.BamReader("hello")

    def test_BamReader_expects_an_index_file(self):
        with self.assertRaises(IOError):
            parser.BamReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/foo.bam"))
            
    def test_fetch_expects_str(self):
        """
        ngsci.reads.BamReader.fetch throws a TypeError when the chrom isn't a str
        """
        with self.assertRaises(TypeError):
            self.bamreader.fetch(None, 0, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(1, 0, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(1.0, 0, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch([], 0, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch({}, 0, 75)

    def test_fetch_expects_known_chrom(self):
        """
        ngsci.reads.BamReader.fetch throws a TypeError when the chrom isn't in the index
        """
        with self.assertRaises(TypeError):
            self.bamreader.fetch("hello", 0, 75)

    def test_fetch_expects_positive_int(self):
        """
        ngsci.reads.BamReader.fetch throws a TypeError when the start site isn't a positive int
        """
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, None, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, "hello", 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 1.0, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, [], 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, {}, 75)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, -1, 75)

    def test_fetch_expects_another_positive_int(self):
        """
        ngsci.reads.BamReader.fetch throws a TypeError when the step size isn't a positive int
        """
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, None)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, "hello")
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 1.0)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, [])
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, {})
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, -1)

    def test_fetch_expects_strand_specific_to_be_a_bool(self):
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific=None)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific="hello")
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific=1)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific=1.0)
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific=[])
        with self.assertRaises(TypeError):
            self.bamreader.fetch(self.chrom, 0, 76, strand_specific={})
            
    def test_get_complexity_returns_tuple_of_ints(self):
        iter = self.samfile.fetch(self.chrom, 0, 75)
        reads = [x for x in iter]
        self.assertEqual(len(reads), 100)

        result = self.bamreader._get_complexity(reads, False)
        self.assertIsInstance(result, tuple) # Returns a tuple
        self.assertEqual(len(result), 2) # of length 2
        self.assertIsInstance(result[0], int) # containing two ints
        self.assertIsInstance(result[1], int)
        self.assertEqual(result[0], 95)
        self.assertEqual(result[1], 0)
            
if __name__ == '__main__':
    unittest.main()
