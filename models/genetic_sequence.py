@dataclass(frozen=True)
class GeneticSequence:
  sequence: str
  is_dna: bool  # Flag to distinguish DNA (True) or RNAa (False)
