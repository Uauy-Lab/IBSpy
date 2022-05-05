from configparser import Error
import logging
import unittest

from IBSpy.IBSpy_options import IBSpyOptions, parse_str_to_list

class TestOptions(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.result = ["cadenza", "robigus"]
        self.options = IBSpyOptions()
        #self.options.metadata = self.samples_metadata
        self.options.chromosomes = "./tests/data/affinity/chromosome_lengths.tsv"

    def test_parse_list(self):
        test_1 = parse_str_to_list(self.result)
        self.assertListEqual(test_1, self.result)
        test_2 = parse_str_to_list("cadenza,robigus")
        self.assertListEqual(test_2, self.result)
        test_3 = parse_str_to_list("tests/data/lines.txt")
        self.assertListEqual(test_3, self.result)
    
    def test_regions(self):
        self.options.region = "chr1A__chi:100-200"
        chr, start, end, asm = self.options.region
        self.assertEqual("chr1A__chi", chr)
        self.assertEqual(100, start)
        self.assertEqual(200, end)
        self.assertEqual("chinese", asm)

        self.options.region = "chr1A__jag:100-200"
        chr, start, end, asm = self.options.region
        self.assertEqual("chr1A__jag", chr)
        self.assertEqual(100, start)
        self.assertEqual(200, end)
        self.assertEqual("jagger", asm)
        self.options.region = "chr1A__jagg:100-200"
        def f():
            self.options.region
        self.assertRaises(ValueError, f)
    
        
if __name__ == '__main__':
    unittest.main()