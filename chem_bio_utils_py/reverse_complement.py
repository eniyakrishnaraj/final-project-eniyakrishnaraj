class ReverseComplement:
  def __init__(self):
    self.valid_nucleotides = set("ATCG")
    self.complement_table = {
      "A": "T",
      "T": "A",
      "C": "G",
      "G": "C"
    }

  def reverseComp(self, seq):
    # Ensures sequence can be handled even if it is in lower case
    sequence = seq.upper()

    # Check if sequence is empty
    if not sequence:
      return "DNA sequence is empty."

    # Check if the sequence contains only valid nucleotides
    if not all(nucleotide in self.valid_nucleotides for nucleotide in sequence):
      return "Invalid nucleotide found in DNA sequence."
      
    complement_seq = ""

    for nucleotide in reversed(sequence):
      complement_seq += self.complement_table.get(nucleotide, "")

    return complement_seq
