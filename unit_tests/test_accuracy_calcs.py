import unittest
import accuracy

class TestAccuracy(unittest.TestCase):

    def test_accuracy_all_good(self):

        called_alleles=[2,3,5]
        expected_alleles=[2,3,5]
        result=accuracy.allele_accuracy(called_alleles, expected_alleles)
        self.assertEqual(result, 100)

        called_alleles=[2,3,5]
        expected_alleles=[2,-1,5]
        result=accuracy.allele_accuracy(called_alleles, expected_alleles)
        self.assertEqual(result, 100)

    def test_accuracy_some_diff(self):

        called_alleles=[2,8,90]
        expected_alleles=[2,10,100]
        result=accuracy.allele_accuracy(called_alleles, expected_alleles)
        self.assertEqual(result,70)

        called_alleles=[2,8, 100]
        expected_alleles=[2,50]
        result=accuracy.allele_accuracy(called_alleles, expected_alleles)
        self.assertAlmostEqual(result, 66.6666, 2)

        called_alleles=[2,8]
        expected_alleles=[2,8,100]
        result=accuracy.allele_accuracy(called_alleles, expected_alleles)
        self.assertAlmostEqual(result, 66.6666, 2)

if __name__ == '__main__':
    unittest.main()
