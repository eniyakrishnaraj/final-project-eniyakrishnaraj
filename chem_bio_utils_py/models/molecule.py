from dataclasses import dataclass

@dataclass(frozen=True)
class MoleculeModel:
  inchi: str
  smiles: str
  molecular_weight: float
  molecular_formula: str
  fingerprint: object
