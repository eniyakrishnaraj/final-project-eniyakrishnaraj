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
    return "Invalid input. Not valid SMILES or InChI."

  # Set the valence state explicitly for each atom to suppress the unusual valence warning
  for atom in mol.GetAtoms():
    if atom.GetTotalValence() != atom.GetExplicitValence():
      atom.SetFormalCharge(atom.GetTotalValence() - atom.GetExplicitValence())
    
  # Processes the input and returns a molecule with all the attributes
  molecule = MoleculeModel(
    inchi=Chem.MolToInchi(mol),
    smiles=Chem.MolToSmiles(mol),
    molecular_weight=Descriptors.MolWt(mol),
    molecular_formula=molecularFormula(mol),
    fingerprint=Chem.RDKFingerprint(mol))
  return molecule

# Returns formula of a given molecule
def molecularFormula(molecule):
  if not molecule:
    return "Input is empty."
  elements = [atom.GetSymbol() for atom in molecule.GetAtoms()]
  molecular_formula = ''.join(sorted(elements))
  return molecular_formula
