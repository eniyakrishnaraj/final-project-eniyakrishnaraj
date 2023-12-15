from chem_bio_utils_py.operations.process_molecule import processMolecule

def process_mol_file(file_path):
  molecule_list = []
  with open(file_path, 'r') as file:
    list = file.read().splitlinese()
    for mol in list:
      molecule = processMolecule(mol)  # Using processMolecule function to process each input
      if isinstance(molecule, str):
        print(molecule)  # Print error message if processing fails
      else:
        molecule_list.append(molecule)  # Add MoleculeModel instance to list
  return molecule_list

if __name__ == "__main__":
  smiles_file_path = 'smiles.txt'
  smiles_molecules = process_mol_file(smiles_file_path)
