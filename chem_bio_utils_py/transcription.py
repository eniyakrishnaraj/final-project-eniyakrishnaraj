class DNATranscriber:
  def __init__(self):
    self.transcription_table = {
      "A": "U",
      "T": "A",
      "C": "G",
      "G": "C"
    }

  def transcribe(self, dna_seq):
    # Ensures sequence can be handled even if it is in lower case
    sequence = dna_seq.upper()

    # Check if sequence is empty
    if not dna_seq:
      return "DNA sequence is empty."

    # Check if the sequence contains only valid nucleotides
    valid_nucleotides = set("ATCG")
    if not all(nucleotide in valid_nucleotides for nucleotide in dnaseq):
      return "Invalid nucleotide found in DNA sequence."
      
    transcribed_sequence = ""

    for nucleotide in sequence:
      transcribed_sequence += self.transcription_table.get(nucleotide, "")

    return transcribed_sequence
