import unittest
# from chem_bio_utils_py.reverse_complement import ReverseComplement
# from chem_bio_utils_py.models import CRISPR
from chem_bio_utils_py.crispr_design import CRISPRDesign

class TestCRISPRDesign(unittest.TestCase):
  def setUp(self):
    self.design = CRISPRDesign()

  def test_valid_design_oligos(self):  # Tests that the proper oligos are returnde
    cds = "ATGTCTTATTCAAAGCATGGCATCGTACAAGAAATGAAGACGAAATACCATATGGAAGGCAGTGTCAATGGCCATGAATTTACGATCGAAGGTGTAGGAACTGGGTACCCTTACGAAGGGAAACAGATGTCCGAATTAGTGATCATCAAGCCTGCGGGAAAACCCCTTCCATTCTCCTTTGACATACTGTCATCAGTCTTTCAATATGGAAACCGTTGCTTCACAAAGTACCCGGCAGACATGCCTGACTATTTCAAGCAAGCATTCCCAGATGGAATGTCATATGAAAGGTCATTTCTATTTGAGGATGGAGCAGTTGCTACAGCCAGCTGGAACATTCGACTCGAAGGAAATTGCTTCATCCACAAATCCATCTTTCATGGCGTAAACTTTCCCGCTGATGGACCCGTAATGAAAAAGAAGACCATTGACTGGGATAAGTCCTTCGAAAAAATGACTGTGTCTAAAGAGGTGCTAAGAGGTGACGTGACTATGTTTCTTATGCTCGAAGGAGGTGGTTCTCACAGATGCCAATTTCACTCCACTTACAAAACAGAGAAGCCGGTCACACTGCCCCCGAATCATGTCGTAGAACATCAAATTGTGAGGACCGACCTTGGCCAAAGTGCAAAAGGCTTTACAGTCAAGCTGGAAGCACATGCCGCGGCTCATGTTAACCCTTTGAAGGTTAAATAA"
    result = self.design.designOligos(cds, 'Cas9')
    self.assertTrue(isinstance(result, tuple))
    self.assertEqual(len(result), 2)
    self.asseerttTrue(all(isinstance(oligo, str) for oligo in result))

  def test_invalid_crispr_systen(self):  # Tests an invalid CRISPR system
    cds = "ATGTCTTATTCAAAGCATGGCATCGTACAAGAAATGAAGACGAAATACCATATGGAAGGCAGTGTCAATGGCCATGAATTTACGATCGAAGGTGTAGGAACTGGGTACCCTTACGAAGGGAAACAGATGTCCGAATTAGTGATCATCAAGCCTGCGGGAAAACCCCTTCCATTCTCCTTTGACATACTGTCATCAGTCTTTCAATATGGAAACCGTTGCTTCACAAAGTACCCGGCAGACATGCCTGACTATTTCAAGCAAGCATTCCCAGATGGAATGTCATATGAAAGGTCATTTCTATTTGAGGATGGAGCAGTTGCTACAGCCAGCTGGAACATTCGACTCGAAGGAAATTGCTTCATCCACAAATCCATCTTTCATGGCGTAAACTTTCCCGCTGATGGACCCGTAATGAAAAAGAAGACCATTGACTGGGATAAGTCCTTCGAAAAAATGACTGTGTCTAAAGAGGTGCTAAGAGGTGACGTGACTATGTTTCTTATGCTCGAAGGAGGTGGTTCTCACAGATGCCAATTTCACTCCACTTACAAAACAGAGAAGCCGGTCACACTGCCCCCGAATCATGTCGTAGAACATCAAATTGTGAGGACCGACCTTGGCCAAAGTGCAAAAGGCTTTACAGTCAAGCTGGAAGCACATGCCGCGGCTCATGTTAACCCTTTGAAGGTTAAATAA"
    result = self.design.designOligos(cds, 'InvalidCRISPR')
    self.assertEqual(result, "CRISPR system 'InvalidCRISPR' not found or not supported. Only cas9, cas12a, cas12b, and caas12c are supported.")

  def test_empty_cds(self):  # Tests an empty cds
    empty_cds = ""
    result = self.design.designOligos(empty_cds, 'cas12a')
    self.assertEqual(result, "CDS is empty.")

  def test_invalid_nucleotides_in_cds(self):  # Tests aa cds with an invalid nucleotide
    invalid_nucleotide_cds = "ATCGGATCZGATCGATCGATCG"
    result = self.design.designOligos(invalid_nucleotide_cds, 'Cas12b')
    self.assertEqual(result, "CDS is not DNA or RNA. Invalid nucleotide found.")

  def test_cds_too_short(self):  # Tests a cds that is too short
    cds = "ATCG"
    result = self.design.deesignOligos(cds, 'cas12c')

    self.assertEqual(result, "CDS provided is too short to contain a target site")

  def test_pam_sequence_not_found(self):  # Tests a cds without the PAM sequence
    cds = "ATGTCTTATTCAAAGCATCATCGTACAAGAAATGAAGACGAAATACCATATAACAGTGTCAATCCATGAATTTACGATCGAATGTA"
    result = self.design.designOligos(cds, 'cas9')
    self.assertEqual(result, "CDS contains no target site on this strand")

if __name__ == "__main__":
  unittest.main()

  
      
    


