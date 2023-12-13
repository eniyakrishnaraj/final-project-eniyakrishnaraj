import unittest
from chem_bio_util_py.translate import Translate

class TestTranslate(unittest.TestCase):
  def setUp(self):
    self.translator = Translate()

  def test_translation(self):  # Test a valid DNA sequence
    dna_sequence = "ATGTCATCGTAG"
    expected_protein_sequence = "MSS*"
    result = self.translator.run(dna_sequence)
    self.assertEqual(result, expected_protein_sequence)

  def test_invalid_sequence(self):
    # Test an invalid DNA sequence (not divisible by 3)
    invalid_dna_sequence = "ATGT"
    result = self.translator.run(invalid_dna_sequence)
    self.assertEqual(result, "")

    # Test an empty DNA sequence
    empty_dna_sequence = ""
    result_empty = self.translator.run(empty_dna_sequence)
    self.assertEqual(result_empty, "")

  def test_invalid_nucleotide(self):  # Test a DNAA sequence with an invalid nucleotide
    invalid_nucleotide_sequence = "ATGTXATCGTAG"
    result_invalid_nucleotide = self.translator.run(invalid_nucleotide_sequence)
    self.assertEqual(result_invalid_nucleotide, "")

if __name__ == "__main__":
  unittest.main()
