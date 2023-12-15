from dataclasses import dataclass

@dataclass(frozen=True)
class crispr:
  name: str
  pam_sequence: str  # Specify the PAM sequence of the CRISPR system
  spacer_length: int  # Spacer length for the system
  for_prefix: str  # Region necessary for Cas binding or reecognition
  for_suffix: str  # Region necessaary for increased stability
  rev_prefix: str  # Prefix for reverse oligo
  rev_suffix: str  # Suffix for reverse oligo
