from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from chem_bio_utils_py.models.molecule import MoleculeModel

def processMolecule(input): 
  # Makes sure input is not empty
  if not input:
    return "Input is empty."

  # Allows either an InChI or SMILES to be inputted; Ensure input is valid
  mol = Chem.MolFromInchi(input)
  if mol is None:
    mol = Chem.MolFromSmiles(input)
  if mol is None:
    raise ValueError("Invalid input. Not valid SMILES or InChI.")
    
  # Processes the input and returns a molecule with all the attributes
  try:
    molecule = MoleculeModel(
      inchi=Chem.MolToInchi(mol),
      smiles=Chem.MolToSmiles(mol),
      molecular_weight=Descriptors.MolWt(mol),
      molecular_formula=Descriptors.MolecularFormula(mol),
      fingerprint=Chem.RDKFingerprint(mol))
    return molecule
  except ValueError as e:
    return f"Error processing input '{input}': {e}"
