from dataclasses import dataclass, field
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

@dataclass(frozen=True)
class MoleculeModel:
  input: str
  inchi: str = field(init=False)  # String notation used as ideentifier for chemical substances
  smiles: str = field(init=False)  # String notation representing the structure of a molecule
  molecular_weight: float = field(init=False)  # Sum of atomic weights of all the atoms in a molecule
  molecular_formula: str = field(init=False)  # Molecular formula indicating types and numbeers of atoms in a molecule
  fingerprint: object = field(init=False)  # Binary representation of molecular structure

  def __post_init__(self):
    mol = Chem.MolFromInchi(self.input)
    if mol is None:
      mol = Chem.MolFromSmiles(self.input)
    if mol is None:
      raise ValueError("Invalid input. Not SMILES or InChI.")

    self.inchi = Chem.MolToInchi(mol)
    self.smiles = Chem.MolToSmiles(mol)
    self.molecular_weight = Descriptors.MolWt(mol)
    self.molecular_formula = Descriptors.MolecularFormula(mol)
    self.fingerprint = Chem.RDKFingerprint(mol)
