"""
Helper functions for design automation and analysis.
"""


from pathlib import Path
from typing import NamedTuple

from biotite.structure.io import load_structure


class Redesign(NamedTuple):
    sequence: list[str]
    merged_sequence: str
    unique_chains: int
    score: float
    recovery: float
    source: str


def get_total_residues(pdb: Path) -> int:
    """
    Get the total number of residues in the structure.
    """
    structure = load_structure(pdb)
    return len(structure[structure.atom_name == "CA"].res_id)


def parse_redesign(annotation: str, sequence: str) -> Redesign:
    """
    From the ProteinMPNN output lines, build the Redesign data.
    """
    sequences = sequence.strip().split("/")

    unique_sequences = []
    for seq in sequences:
        if seq not in unique_sequences:
            unique_sequences.append(seq)
    merged_sequence = "_".join(unique_sequences)

    unique_chains = len(set(sequences))

    score = [
        float(info.split("=")[1])
        for info in annotation.split(",")
        if "score" in info.split("=")[0]
    ][0]

    if "seq_recovery" in annotation:
        recovery = [
            float(info.split("=")[1])
            for info in annotation.split(",")
            if "seq_recovery" in info.split("=")[0]
        ][0]
    else:
        recovery = 1.0

    if "sample" in annotation:
        source = "design"
    else:
        source = "wt"

    return Redesign(
        sequence=sequences,
        merged_sequence=merged_sequence,
        unique_chains=unique_chains,
        score=score,
        recovery=recovery,
        source=source,
    )


def load_proteinmpnn_output(fasta: Path) -> list[Redesign]:
    """
    Read a ProteinMPNN output file and return a list of designs.
    """
    with open(fasta, mode="r", encoding="utf-8") as fa_file:
        lines = fa_file.readlines()

    redesigns = []
    for i in range(0, len(lines), 2):
        redesign = parse_redesign(lines[i], lines[i + 1])
        if redesign not in redesigns:
            redesigns.append(parse_redesign(lines[i], lines[i + 1]))

    return redesigns
