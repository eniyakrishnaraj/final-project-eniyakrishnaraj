from chem_bio_utils_py.models.molecule import MoleculeModel

def processMolecule(input):  # Can take in either SMILES or InChI
  # Check that input is not empty string
  if not input:
    return "Input is empty."

  # Processes the input and returns a molecule with all the attributes
  try:
    molecule = MoleculeModel(input)
    return molecule
  except ValueError as e:
    return f"Error processing input '{input}': {e}"
