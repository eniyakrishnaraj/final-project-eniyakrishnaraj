import unittest
from chem_bio_utils_py.operations.process_molecules import processMolecule

class TestProcessMolecule(unittest.TestCase):
  def setUp(self):
    smiles = "[H-]"
    inchi = "InChI=1S/H/q-1"
    
  def test_valid_smiles(self):  # Test a valid SMILES string
    smiles = "[H-]"
    inchi = "InChI=1S/H/q-1"
    result = processMolecule(smiles)
    self.assertIsNotNone(result)
    self.assertNotIsInstance(result, str)
    self.assertEqual(inchi, result.inchi)

  def test_invalid_smiles(self):  # Test an invalid SMILES string
    invalid_smiles = "InvalidSMILES"
    result = processMolecule(invalid_smiles)
    self.assertIsInstance(result, str)
    self.assertEqual(result, "Error processing input 'InvalidSMILES': Invalid input. Not SMILES or InChI.")

  def test_valid_inchi(self):  # Test a valid InChI string
    inchi = "InChI=1S/H/q-1"
    smiles = "[H-]"
    result = processMolecule(inchi)
    self.assertIsNotNone(result)
    self.assertNotIsInstance(result, str)
    self.assertEqual(smiles, result.smiles)

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
