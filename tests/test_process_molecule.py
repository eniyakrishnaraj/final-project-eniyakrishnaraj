import unittest
from chem_bio_utils_py.operations.process_molecules import processMolecule

class TestProcessMolecule(unittest.TestCase):
  def test_valid_smiles(self):  # Test a valid SMILES string
    input = "[CH4+]"
    result_inchi = "InChI=1S/CH4/h1H4/q+1"
    result_mw = 16.043
    result_mf = "CH4"
    result = processMolecule(input)
    self.assertEqual(result_inchi, result.inchi)
    self.assertEqual(result_mw, result.molecular_weight)
    self.assertEqual(result_mf, result.molecular_formula)

  def test_invalid_smiles(self):  # Test an invalid SMILES string
    invalid_smiles = "InvalidSMILES"
    result = processMolecule(invalid_smiles)
    self.assertEqual(result, "Invalid input. Not valid SMILES or InChI.")

  def test_valid_inchi(self):  # Test a valid InChI string
    input = "InChI=1S/CH4/h1H4/q+1"
    result_smiles = "[CH4+]"
    result_mw = 16.043
    result_mf = "CH4"
    result = processMolecule(input)
    self.assertEqual(result_smiles, result.smiles)
    self.assertEqual(result_mw, result.molecular_weight)
    self.assertEqual(result_mf, result.molecular_formula)

  def test_invalid_inchi(self):  # Test an invalid InChI string
    invalid_inchi = "InvalidInChI"
    result = processMolecule(invalid_inchi)
    self.assertEqual(result, "Invalid input. Not valid SMILES or InChI.")

  def test_empty_input(self):  # Test an empty string
    empty_input = ""
    result = processMolecule(empty_input)
    self.assertIsInstance(result, str)
    self.assertEqual(result, "Input is empty.")

if __name__ == '__main__':
  unittest.main()
