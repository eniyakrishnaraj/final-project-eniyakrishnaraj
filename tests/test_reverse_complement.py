import unittest
from chem_bio_utils_py.operations.reverse_complement import ReverseComplement

class TestReverseComplement(unittest.TestCase):
  def setUp(self):
    self.revcomp = ReverseComplement()

  def test_transcribe(self):  # Tests a valid DNA seequence
    dna_sequence = "ATGTCATTCGTAG"
    expected_comp_sequence = "CTACGAATGACAT"
    result = self.revcomp.reverseComp(dna_sequence)
    self.assertEqual(result, expected_comp_sequence)

  def test_lowercase_sequence(self):  # Test a DNA sequence with lowercase letters
    lowercase_dna_sequence = "atgtcattcgtag"
    expected_comp_sequence = "CTACGAATGACAT"
    result_lowercase = self.revcomp.reverseComp(lowercase_dna_sequence)
    self.assertEqual(result_lowercase, expected_comp_sequence)
    
  def test_empty_sequence(self):  # Test an empty DNA sequence
    empty_dna_sequence = ""
    result_empty =  self.revcomp.reverseComp(empty_dna_sequence)
    self.assertEqual(result_empty, "DNA sequence is empty.")

  def test_invalid_nucleotide(self):   # Tests a DNA sequence with invalid nucleotides (not 'A', 'T', 'C', or 'G')
    invalid_nucleotide_sequence = "ATGTXATCGTAG"
    result_invalid_nucleotide = self.revcomp.reverseComp(invalid_nucleotide_sequence)
    self.assertEqual(result_invalid_nucleotide, "Invalid nucleotide found in DNA sequence.")

if __name__ == "__main__":
  unittest.main()
