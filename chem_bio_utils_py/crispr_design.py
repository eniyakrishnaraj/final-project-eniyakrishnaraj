from chem_bio_utils_py.models import CRISPR

class CRISPRDesign:
  def __init__(self, target_sequence, organism_genome):
    self.target_sequence = target_sequence
    self.organism_genome = organism_genome
    self.valid_crispr_systems = self.initialize.crispr_systems()

  def initialize_crispr_systems(self):
    cas9 = CRISPR(name='Cas9', pam_sequence='GG')  # NGG where N can be any nucleotide
