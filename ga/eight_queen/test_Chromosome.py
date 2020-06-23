import unittest

from Chromosome import Chromosome


class ChromosomeTestCase(unittest.TestCase):
    def test_constructor_invalid_element_type(self):
        try:
            Chromosome(['a', 'b', 'c', 'd'])
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_constructor_negative_integer(self):
        try:
            Chromosome([-1, -2, 4, 5, 6])
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_constructor_valid_data(self):
        try:
            Chromosome([0, 1, 2, 3, 4, 5, 6])
            pass
        except TypeError:
            self.fail("TypeError exception raised with valid data.")

    def test_constructor_valid_data_invalid_length(self):
        try:
            Chromosome([0])
        except ValueError:
            pass
        else:
            self.fail("Unknown Exception raised.")


if __name__ == '__main__':
    unittest.main()
