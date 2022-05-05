import logging
import unittest

from IBSpy.IBSpy_options import parse_str_to_list

class TestOptions(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.result = ["cadenza", "robigus"]

    def test_parse_list(self):
        test_1 = parse_str_to_list(self.result)
        self.assertListEqual(test_1, self.result)
        test_2 = parse_str_to_list("cadenza,robigus")
        self.assertListEqual(test_2, self.result)
        test_3 = parse_str_to_list("tests/data/lines.txt")
        self.assertListEqual(test_3, self.result)
    
    def test_regions(self):
        pass 

if __name__ == '__main__':
    unittest.main()