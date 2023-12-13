class Translate:
  def __init__(self):
    self.codon_table = {
      "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
      "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
      "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
      "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
      "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
      "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
      "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
      "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
      "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
      "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
      "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
      "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
      "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
      "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
      "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
      "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G"
    }

  def run(self, seq):
    dnaseq = seq.upper()

    # Check if sequence is empty or not divisible by 3 (to be split into codons)
    if not dnaseq or len(dnaseq) % 3 != 0:
      print("Dna sequence is not valid.")

    # Initializes empty string for protein sequence
    proteinseq = ""
    # Populates protein sequence with appropriate aamino acid from codon_table
    for i in range(0, len(dnaseq), 3):
      codon = dnaseq[i:i + 3]
      proteinseq += self.codon_table.get(codon, "")

    return proteinseq
