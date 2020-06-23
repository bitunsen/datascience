import unittest

from Gene import Gene


class GeneTestCase(unittest.TestCase):

    def test_constructor(self):
        try:
            gene_1 = Gene(4, "5")
        except TypeError:
            pass
        else:
            self.fail("Unknown Exception raised.")

    def test_equals_true(self):
        gene_1 = Gene(4, 5)
        gene_2 = Gene(4, 5)
        expected = gene_1.equals(gene_2)
        self.assertEqual(expected, True)

    def test_equals_false(self):
        gene_1 = Gene(4, 5)
        gene_2 = Gene(4, 8)
        expected = gene_1.equals(gene_2)
        self.assertEqual(expected, False)

    def test_equals_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.equals, "String")

    def test_prints_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.print, 123)

    def test_prints_type_check(self):
        gene_1 = Gene(4, 5)
        self.assertRaises(TypeError, gene_1.print, 123)


if __name__ == '__main__':
    unittest.main()
