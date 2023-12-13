from dataclasses import dataclass

@dataclass(frozen=True)
class crispr:
  name: str
  pam_sequence: str
