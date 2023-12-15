from chem_bio_utils_py.models.crispr import crispr
from chem_bio_utils_py.reverse_complement import ReverseComplement

class CRISPRDesign:
  def __init__(self):
    self.valid_crispr_systems = self.initialize_crispr_systems()
    self.valid_nucleotides = set("ATCGU")
    self.revcomp = ReverseComplement()

  """This method is meant to initialize the CRISPR systtems. All the attributes were designed based
     on literature review (sources below). This initialization can obviously be added upon as 
     more research is conducted on more CRISPR systems. https://star-protocols.cell.com/protocols/1707"""
  def initialize_crispr_systems(self):
    # PAM: NGG, GG for simplicity; Simple, efficient and open-source CRISPR/Cas9 strategy for multi-site genome editing in Populus tremula × alba: https://academic.oup.com/treephys/article/41/11/2216/6270869
    cas9 = crispr(name='Cas9', pam_sequence='GG', spacer_length=20, for_prefix='AATGGTCTCAA', for_suffix='GCTTAGAGACCAAT', rev_prefix='CTCACTAG', rev_suffix='TAG')
    # PAM: TTTV, TTT for simplicity; Repeat expansion and methylation state analysis with nanopore sequencing: https://www.biorxiv.org/content/10.1101/480285v1.full
    cas12a = crispr(name='Cas12a', pam_sequence='TTT', spacer_length=30, for_prefix='CGGCAGCC', for_suffix='GGAAAAAC', rev_prefix='TTAGCGGCTT', rev_suffix='GGCAAAAGAG')
    # PAM: TCCC; Engineering of CRISPR-Cas12b for human genome editing: https://www.nature.com/articles/s41467-018-08224-4#MOESM1
    cas12b = crispr(name='Cas12b', pam_sequence='TCCC', spacer_length=30, for_prefix='ATTNT', for_suffix='CTT', rev_prefix='GAGCTGCG', rev_suffix='CTGTGATAC')
    # PAM: TTN, TT for simplicity; A naturally DNase-free CRISPR-Cas12c enzyme silences gene expression: https://www.sciencedirect.com/science/article/pii/S1097276522003835#mmc3
    cas12c = crispr(name='Cas12c', pam_sequence='TT', spacer_length=30, for_prefix='GCCTGCCCGCAGA', for_suffix='CTTTCCAGTC', rev_prefix='TGGCTGT', rev_suffix='ACCGTAA')

    return [cas9, cas12a, cas12b, cas12c]

  def designOligos(self, cds, crispr_system_name):
    # Find the CRISPR system by name
    crispr = next((c for c in self.valid_crispr_systems if c.name.lower() == crispr_system_name.lower()), None)
    if not crispr:
      return f"CRISPR system '{crispr_system_name}' not found or not supported. Only cas9, cas12a, cas12b, and caas12c are supported."
      
    # Check that cds is not an empty string
    if not cds:
      return "CDS is empty."

    # Check if the length of the cds is long enough based on the spacer length and PAM sequence
    if len(cds) < (crispr.spacer_length + len(crispr.pam_sequence)):
      return "CDS provided is too short to contain a target site."

    # Check that the cds only has appropriate nucleotides
    cds = cds.upper()
    if not all(nucleotide in self.valid_nucleotides for nucleotide in cds):
      return "CDS is not DNA or RNA. Invalid nucleotide found."

    # Calculate the start index for the PAM sequence
    PAM_start = cds.find(crispr.pam_sequence, crispr.spacer_length + len(crispr.pam_sequence)) - 1
    if PAM_start < 0:
      return "CDS contains no target site on this strand"

    # Construct the spacer
    spacer = cds[PAM_start - crispr.spacer_length:PAM_start]
    revCompSpacer = self.revcomp.reverseComp(spacer)

    for_oligo = crispr.for_prefix + spacer + crispr.for_suffix
    rev_oligo = crispr.rev_prefix + revCompSpacer + crispr.rev_suffix

    return [for_oligo, rev_oligo]
