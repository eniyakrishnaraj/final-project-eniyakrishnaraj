import unittest
from chem_bio_utils_py.transcription import DNATranscriber

class TestDNATranscriber(unittest.TestCase):
  def setUp(self):
    self.transcriber = DNATranscriber()

  def test_transcribe(self):  # Tests a valid DNA seequence
    dna_sequence = "ATGTCATTCGTAG"
    expected_rna_sequence = "UACAGUAAGCAUC"
    result = self.transcriber.transcribe(dna_sequence)
    self.assertEqual(result, expected_rna_sequence)

  def test_lowercase_sequence(self):  # Test a DNA sequence with lowercase letters
    lowercase_dna_sequence = "atgtcattcgtag"
    expected_rna_sequence = "UACAGUAAGCAUC"
    result_lowercase = self.transcriber.transcribe(lowercase_dna_sequence)
    self.assertEqual(result_lowercase, expected_rna_sequence)
    
  def test_empty_sequence(self):  # Test an empty DNA sequence
    empty_dna_sequence = ""
    result_empty =  self.transcriber.transcribe(empty_dna_sequence)
    self.assertEqual(result_empty, "DNA sequence is empty.")

  def test_invalid_nucleotide(self):   # Tests a DNA sequence with invalid nucleotides (not 'A', 'T', 'C', or 'G')
    invalid_nucleotide_sequence = "ATGTXATCGTAG"
    result_invalid_nucleotide = self.transcriber.transcribe(invalid_nucleotide_sequence)
    self.assertEqual(result_invalid_nucleotide, "Invalid nucleotide found in DNA sequence.")

if __name__ == "__main__":
  unittest.main()
