### Eniya Krishnaraj BioE134 Individual Project
# ChemBioUtilsPy
ChemBioUtilsPy is a Python project designed to have some chemical and biological operations that are useful. It has some basic operations for translation, transcription, and reverse complementation. Then, it has a Molecule Model that has a SMILES, InChI, molecular weight, and molecular formula. This Model is utilized in a function that takes in input (either SMILES or InChI), and returns a MoleculeModel with all the attributes. This function was then implemented in order to take in a txt file with a list of SMILES, process them, and return a list of instances of MoleculeModel. Finally, this project designs oligos for 4 CRISPR-Cas systems: Cas9, Cas12a, Cas12b, and Cas12c. It utilizes the crispr Model, and can easily be expanded on with more Cas systems. The oligos' "templates" (prefix and suffix) came from extensive literature review.

This project follows the Function/Model pattern, separating data models from business logic. It utilizes unittest as a testing framework for all the functions.

# Author
Eniya Krishnaraj - EK - eniyakrishnaraj

# Acknowledgements
Thank you to Prof Anderson and Riddhi for teaching both the biological and coding knowledge needed. 
