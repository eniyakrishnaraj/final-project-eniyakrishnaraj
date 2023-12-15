from dataclasses import dataclass

@dataclass(frozen=True)
class GeneticSequence:
  sequence: str
  is_dna: bool # Flag to distinguish DNAA (True) or RNA (False)
