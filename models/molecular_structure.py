@dataclass(frozen=True)
class MolecularStructure:
  smiles: str
  inchi: str
  molecular_formula: str
  molecular_weight: float
  charge: int
  stereochemistry: str
  bond_info: dict
  structuraal_isomers: list # List of structural isomers
