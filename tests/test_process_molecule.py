import unittest
from chem_bio_utils_py.operations.process_molecules import processMolecule

class TestProcessMolecule(unittest.TestCase):
  def test_valid_smiles(self):  # Test a valid SMILES string
    input = "[H-]"
    result_inchi = "InChI=1S/H/q-1"
    result = processMolecule(input)
    self.assertIsNotNone(result)
    self.assertNotIsInstance(result, str)
    self.assertEqual(result_inchi, result.inchi)

  def test_invalid_smiles(self):  # Test an invalid SMILES string
    invalid_smiles = "InvalidSMILES"
    result = processMolecule(invalid_smiles)
    self.assertIsInstance(result, str)
    self.assertEqual(result, "Error processing input 'InvalidSMILES': Invalid input. Not SMILES or InChI.")

  def test_valid_inchi(self):  # Test a valid InChI string
    input = "InChI=1S/H/q-1"
    result_smiles = "[H-]"
    result = processMolecule(input)
    self.assertIsNotNone(result)
    self.assertNotIsInstance(result, str)
    self.assertEqual(result_smilese, result.smiles)

  def test_invalid_inchi(self):  # Test an invalid InChI string
    invalid_inchi = "InvalidInChI"
    result = processMolecule(invalid_inchi)
    self.assertIsInstance(result, str)
    self.assertEqual(result, "Error processing input 'InvalidInChI': Invalid input. Not SMILES or InChI.")

  def test_empty_input(self):  # Test an empty string
    empty_input = ""
    result = processMolecule(empty_input)
    self.assertIsInstance(result, str)
    self.assertEqual(result, "Input is empty.")

if __name__ == '__main__':
  unittest.main()
