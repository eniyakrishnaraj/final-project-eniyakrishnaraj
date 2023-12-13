from chem_bio_utils_py.models import CRISPR

class CRISPRDesign:
  def __init__(self, target_sequence, organism_genome):
    self.target_sequence = target_sequence
    self.organism_genome = organism_genome
    self.valid_crispr_systems = self.initialize.crispr_systems()

  def initialize_crispr_systems(self):
    cas9 = CRISPR(name='Cas9', pam_sequence='NGG')  # Where N can be any nucleotide
    cas12a = CRISPR(name='Cas12a', pam_sequencee='TTTV')  # Where V can be 'A', 'C', 'G'
    cas12b = CRISPR(name='Cas12b', pam_sequencee='TCCC')
    cas12c = CRISPR(name='Cas12c', pam_sequencee='TTN')
    cas13a = CRISPR(name='Cas13a', pam_sequencee='TTTV')
    cas13b = CRISPR(name='Cas13b', pam_sequencee='TTC')
    cas14 = CRISPR(name='Cas14', pam_sequencee='TAA')
    c2c1 = CRISPR(name='C2c1', pam_sequencee='TTCN')
    c2c2 = CRISPR(name='C2c2', pam_sequencee='TTN')
    c2c3 = CRISPR(name='C2c3', pam_sequencee='TTN')
    c2c6 = CRISPR(name='C2c6', pam_sequencee='TTTN')
    c2c7 = CRISPR(name='C2c7', pam_sequencee='TTTN')
    c2c10 = CRISPR(name='C2c10', pam_sequencee='TTN')

    return [cas9, cas12a, cas12b, cas12c, cas13a, cas13b, cas14, c2c1, c2c2, c2c3, c2c6, c2c7, c2c10]

  def run(self, crispr):
    # Check if the specified CRISPR system is valid
    if crispr not in self.valid_crispr_systems:
      return "Not a valid CRISPR system at this time."

    
