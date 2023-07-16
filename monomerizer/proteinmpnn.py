"""
Functions for making ProteinMPNN input or analyzing output.
"""


from pathlib import Path
from typing import NamedTuple

import biotite.structure as bts

from monomerizer import pdb


class Redesign(NamedTuple):
    sequence: list[str]
    merged_sequence: str
    unique_chains: int
    score: float
    recovery: float
    source: str


def make_fixed_positions_jsonl(
    pdb_file: Path, structure: bts.AtomArray, generated_resids: list[int]
) -> dict[str, dict[str, list[int]]]:
    """
    Given the original unified structure and that with linkers generated,
    determine which residues are in contact with the generated linker and make a .jsonl
    file for ProteinMPNN that would fix all other residues.

    format: {"PDB.stem": {"ChainID": [residue numbers of all fixed]}}
    """
    variable_resids = (
        pdb.get_resids_in_contact(structure, generated_resids) + generated_resids
    )
    fixed_resids = [
        resid
        for resid in structure[structure.atom_name == "CA"].res_id
        if resid not in variable_resids
    ]

    fixed_positions_dict = {
        pdb_file.stem: {structure.chain_id[0]: [int(res_id) for res_id in fixed_resids]}
    }

    return fixed_positions_dict


def make_symmetric_positions_jsonl(
    pdb_file: Path, structure: bts.AtomArray, generated_resids: list[int]
) -> dict[str, dict[str, list[int]]]:
    """
    Given the original unified structure and that with linkers generated,
    determine which residues are in contact with the generated linker and make a .jsonl
    file for ProteinMPNN that would fix all other residues.

    format: {"PDB.stem": {"ChainID": [residue numbers tied position 1]}, {"ChainID: [residue numbers tied position 2]}}
    """
    symmetric_positions_dict = {pdb_file.stem: []}

    resids = structure[structure.atom_name == "CA"].res_id
    original_resids = [resid for resid in resids if resid not in generated_resids]

    original_symmetric_units = []
    last_resid = None
    current_unit = []
    for resid in original_resids:
        if (last_resid is not None) and (resid - 1 > last_resid):
            original_symmetric_units.append(current_unit)
            current_unit = []
        current_unit.append(int(resid))
        last_resid = resid
    original_symmetric_units.append(current_unit)

    generated_symmetric_units = []
    last_resid = None
    current_unit = []
    for resid in generated_resids:
        if (last_resid is not None) and (resid - 1 > last_resid):
            generated_symmetric_units.append(current_unit)
            current_unit = []
        current_unit.append(int(resid))
        last_resid = resid
    generated_symmetric_units.append(current_unit)
    if len(generated_symmetric_units) < len(original_symmetric_units):
        generated_symmetric_units.append([])

    full_symmetric_units = [
        original_symmetric_units[idx] + generated_symmetric_units[idx]
        for idx in range(len(original_symmetric_units))
    ]

    for idx in range(len(full_symmetric_units[0])):
        symmetric_positions = [
            symmetric_unit[idx]
            for symmetric_unit in full_symmetric_units
            if idx < len(symmetric_unit)
        ]
        symmetric_positions_dict[pdb_file.stem].append(
            {structure.chain_id[0]: symmetric_positions}
        )

    return symmetric_positions_dict


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
