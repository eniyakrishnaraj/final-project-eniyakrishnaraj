class DNATranscriber:
  def __init__(self):
    self.transcription_table = {
      "A": "U",
      "T": "A",
      "C": "G",
      "G": "C"
    }

  def transcribe(self, dna_seq):
    if not dna_seq:
      return "DNA sequence is empty."

    valid_nucleotides = set("ATCG")
    if not set(dna_seq).issubset(valid_nucleotides):
      return "Invalid nucleotide found in DNA sequence."

    transcribed_sequence = ""

    for nucleotide in dna_seq.upper():
      transcribed_sequence.add(self.transcription_table.get(nucleotide, ""))

  return transcribed_sequence
