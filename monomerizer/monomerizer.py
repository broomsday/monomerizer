"""
Functions to aid in monomerizing a homo-oligomer.
"""


from pathlib import Path
from dataclasses import dataclass
import subprocess

import numpy as np
from tqdm import tqdm
import biotite.structure as bts
from biotite.structure.io import load_structure, save_structure
import jsonlines

from monomerizer.utils import get_total_residues
from monomerizer.pdb import revert_to_input, compute_symmetry_rmsd
from monomerizer.constants import (
    ANGSTROMS_PER_HELICAL_RESIDUE,
    PEPTIDE_BOND_ATOMS,
    UPPER_HELICAL_MODIFIER,
    LOWER_HELICAL_MODIFIER,
    EXTENDED_RESIDUES_PER_LOOP,
    LOOP_DISTANCE,
)


@dataclass
class GeneratedStructure:
    """Information on an RFDiffusion generated structure."""

    pdb: Path
    total_length: int
    linker_length: int
    symmetry_rmsd: float

    def __init__(
        self, wt_pdb: Path, input_pdb: Path, generated_pdb: Path, linker_length: int
    ):
        wt_structure = load_structure(wt_pdb)
        input_structure = load_structure(input_pdb)
        generated_structure = load_structure(generated_pdb)

        input_resids, generated_resids = get_unified_resids(
            input_structure, generated_structure, linker_length
        )

        reverted_structure = revert_to_input(
            wt_structure,
            generated_structure,
            input_resids,
            generated_resids,
        )
        save_structure(generated_pdb, reverted_structure)

        fixed_positions = make_fixed_positions_jsonl(
            generated_pdb, generated_structure, generated_resids
        )
        with open(generated_pdb.with_suffix(".fixed.jsonl"), mode="wb") as fp:
            with jsonlines.Writer(fp) as writer:
                writer.write(fixed_positions)

        symmetric_positions = make_symmetric_positions_jsonl(
            generated_pdb, generated_structure, generated_resids
        )
        with open(generated_pdb.with_suffix(".symmetric.jsonl"), mode="wb") as fp:
            with jsonlines.Writer(fp) as writer:
                writer.write(symmetric_positions)

        self.pdb = generated_pdb
        self.total_length = get_total_residues(generated_pdb)
        self.linker_length = linker_length
        self.symmetry_rmsd = compute_symmetry_rmsd(reverted_structure, generated_resids)


def compute_res_range(n_nitrogen: bts.Atom, c_carbon: bts.Atom) -> tuple[int, int]:
    """
    Given the nitrogen atom on the N-terminus and the carbon atom on the C-terminus, compute the distance
    and then estimate the min and max number of residues to span that distance, assuming largely helical.
    """
    distance = bts.distance(n_nitrogen, c_carbon)

    max_residues = (
        round(distance / ANGSTROMS_PER_HELICAL_RESIDUE) * UPPER_HELICAL_MODIFIER
        + EXTENDED_RESIDUES_PER_LOOP
    )

    min_residues = (
        round((distance - LOOP_DISTANCE) / ANGSTROMS_PER_HELICAL_RESIDUE)
        * LOWER_HELICAL_MODIFIER
        + EXTENDED_RESIDUES_PER_LOOP
    )

    return (
        max(1, min_residues),
        max_residues,
    )


def get_res_range(structure: bts.AtomArray) -> list[int]:
    """
    Determine the average min and max number of residues to span the chain gaps.
    """
    peptide_bond_atoms = structure[np.isin(structure.atom_name, PEPTIDE_BOND_ATOMS)]

    min_res_values = []
    max_res_values = []

    last_atom = peptide_bond_atoms[0]
    for atom in peptide_bond_atoms:
        if (atom.res_id - last_atom.res_id) > 1:
            min_res, max_res = compute_res_range(atom, last_atom)
            min_res_values.append(min_res)
            max_res_values.append(max_res)
        last_atom = atom
    min_res, max_res = compute_res_range(peptide_bond_atoms[0], last_atom)

    return list(range(round(np.mean(min_res_values)), round(np.mean(max_res_values))))


def parse_manual_res_range(manual_lengths: str) -> list[int]:
    """
    Given a user input for a residue length range, return the list of lengths.
    """
    lengths = []
    for element in manual_lengths.split(","):
        if "-" in element:
            for length in range(
                int(element.split("-")[0]), int(element.split("-")[1]) + 1
            ):
                lengths.append(length)
        else:
            lengths.append(int(element))

    return lengths


def generate_contigs_string(structure: bts.AtomArray) -> str:
    """
    Based on the spacing between N- and C-termini of adjacent chains,
    generates a `contigs` string for use in RFDiffusion that should be
    capable of generating a helical segment to span the termini.
    """
    peptide_bond_atoms = structure[np.isin(structure.atom_name, PEPTIDE_BOND_ATOMS)]

    last_atom = peptide_bond_atoms[0]
    contig_string = f"{last_atom.chain_id}{last_atom.res_id}-"
    for atom in peptide_bond_atoms:
        # when we hit a chain break, add our notation
        if (atom.res_id - last_atom.res_id) > 1:
            contig_string += f"{last_atom.res_id}"
            contig_string += "/XXX_LINKER/"
            contig_string += f"{atom.chain_id}{atom.res_id}-"
        last_atom = atom

    contig_string += f"{last_atom.res_id}"
    contig_string += "/XXX_LINKER"

    return contig_string


def get_shell_script_rfdiffusion(
    input_pdb: Path,
    output_dir: Path,
    res_length: int,
    contigs: str,
    num_designs: int,
    steps: int,
) -> str:
    """
    Write a shell script that will run RFDiffusion using the generated output file and computed contigs.
    """
    shell_str = "#!/bin/bash\n\n"
    shell_str += "cd /home/broom/AlphaCarbon/software/RFdiffusion\n"
    shell_str += "source /usr/etc/profile.d/conda.sh\n"
    shell_str += "conda deactivate\n"
    shell_str += "conda activate SE3nv\n"

    contig_map = contigs.replace("XXX_LINKER", str(res_length))

    shell_str += f"./scripts/run_inference.py "
    shell_str += f"'contigmap.contigs=[{contig_map}]' "
    shell_str += f"inference.num_designs={num_designs} "
    shell_str += f"inference.output_prefix={output_dir.absolute()}/{input_pdb.stem}_{res_length}_monomer "
    shell_str += f"inference.input_pdb={input_pdb.absolute()} "
    shell_str += f"diffuser.T={steps}\n"
    shell_str += "cd /home/broom/AlphaCarbon/code/porepep\n"

    return shell_str


def get_shell_script_pgenerator(
    input_pdb: Path,
    output_dir: Path,
    res_length: int,
    contigs: str,
    num_designs: int,
    steps: int,
) -> str:
    """
    Write a shell script that will run RFDiffusion using the generated output file and computed contigs.
    """
    shell_str = "#!/bin/bash\n\n"
    shell_str += "cd /home/broom/AlphaCarbon/software/protein_generator\n"
    shell_str += "source /usr/etc/profile.d/conda.sh\n"
    shell_str += "conda deactivate\n"
    shell_str += "conda activate proteingenerator\n"

    symmetry = contigs.count("XXX_LINKER")
    contig_map = contigs.replace("XXX_LINKER", str(res_length))
    contig_map = contig_map.replace(
        "/", ","
    )  # modify the RFDiffusion default for Protein Generator

    shell_str += f"python ./inference.py "
    shell_str += f"--contigs {contig_map} "
    shell_str += f"--symmetry {symmetry} "
    shell_str += f"--num_designs {num_designs} "
    shell_str += f"--out {output_dir.absolute()}/{input_pdb.stem}_{res_length}_monomer "
    shell_str += "--save_best_plddt "
    shell_str += f"--pdb {input_pdb.absolute()} "
    shell_str += f"--T {steps}\n"
    shell_str += "cd /home/broom/AlphaCarbon/code/porepep\n"

    return shell_str


def format_design_idx(design_idx: int, method: str):
    """
    Format the design index depending on if this is from RFDiffusion or ProteinGenerator.
    """
    if method == "ProteinGenerator":
        return f"{design_idx:06d}"

    return str(design_idx)


def generate_linked_structures(
    wt_pdb: Path,
    input_pdb: Path,
    output_dir: Path,
    res_lengths: list[int],
    num_designs: int,
    contigs: str,
    method: str = "RFDiffusion",
    steps: int = 25,
) -> list[GeneratedStructure]:
    """
    Run RFDiffusion or ProteinGenerator to generate linker residues across a series of lengths.
    """
    scans = []
    for res_length in tqdm(res_lengths, desc="Generating structures"):
        scan_dir = output_dir / str(res_length)
        scan_dir.mkdir(exist_ok=True, parents=True)

        # if we've already run this length, move to the next
        if len(list(scan_dir.glob("*.pdb"))) < num_designs:
            if method == "ProteinGenerator":
                shell_script = get_shell_script_pgenerator(
                    input_pdb, scan_dir, res_length, contigs, num_designs, steps
                )
            else:
                shell_script = get_shell_script_rfdiffusion(
                    input_pdb, scan_dir, res_length, contigs, num_designs, steps
                )

            subprocess.run(shell_script, shell=True)  # TODO: hide stdout

        for design_idx in range(num_designs):
            run_pdb = (
                scan_dir
                / f"{input_pdb.stem}_{res_length}_monomer_{format_design_idx(design_idx, method)}.pdb"
            )
            scans.append(GeneratedStructure(wt_pdb, input_pdb, run_pdb, res_length))

    return scans


def get_unified_resids(
    input_structure: bts.AtomArray,
    generated_structure: bts.AtomArray,
    linker_length: int,
) -> tuple[list[int], list[int]]:
    """
    Given the input structure and RFD generated structure, return a list of the generated resids.
    """
    # RFDiffusion/ProteinGenerator renumber residues to start from 1 and be completely sequential
    #   so adjust resids to match
    input_resids = input_structure[input_structure.atom_name == "CA"].res_id

    unified_resids = []
    last_resid = None
    chain_breaks = 0
    for resid in input_resids:
        if last_resid == None:
            offset_resid = resid
        elif resid - last_resid > 1:
            chain_breaks += 1
            offset_resid = resid - len(unified_resids) - (linker_length * chain_breaks)

        unified_resids.append(resid - offset_resid + 1)
        last_resid = resid

    # get the residue ids that were NOT in the input structure, these are the generated ones
    generated_resids = [
        res_id
        for res_id in generated_structure[generated_structure.atom_name == "CA"].res_id
        if res_id not in unified_resids
    ]

    return unified_resids, generated_resids


def get_resids_in_contact(
    structure: bts.AtomArray, query_resids: list[int], cutoff_distance: float = 6.0
) -> list[int]:
    """
    Given a structure and query list of residue IDs, return all resids with CA within a cutoff distance.
    Note: Does not include the query list itself.
    """
    cell_list = bts.CellList(
        structure[structure.atom_name == "CA"], cell_size=cutoff_distance
    )
    adjacency_matrix = cell_list.create_adjacency_matrix(cutoff_distance)

    query_res_idxs = [
        resid - 1 for resid in query_resids
    ]  # adjacency matrix is 0-indexed, whereas res_ids are 1-indexd

    contacting_resids = []
    for res_id in structure[structure.atom_name == "CA"].res_id:
        query_adjacency = adjacency_matrix[res_id - 1, query_res_idxs]
        if (True in query_adjacency) and (res_id not in query_resids):
            contacting_resids.append(res_id)

    return contacting_resids


def make_fixed_positions_jsonl(
    pdb: Path, structure: bts.AtomArray, generated_resids: list[int]
) -> dict[str, dict[str, list[int]]]:
    """
    Given the original unified structure and that with linkers generated,
    determine which residues are in contact with the generated linker and make a .jsonl
    file for ProteinMPNN that would fix all other residues.

    format: {"PDB.stem": {"ChainID": [residue numbers of all fixed]}}
    """
    variable_resids = (
        get_resids_in_contact(structure, generated_resids) + generated_resids
    )
    fixed_resids = [
        resid
        for resid in structure[structure.atom_name == "CA"].res_id
        if resid not in variable_resids
    ]

    fixed_positions_dict = {
        pdb.stem: {structure.chain_id[0]: [int(res_id) for res_id in fixed_resids]}
    }

    return fixed_positions_dict


def make_symmetric_positions_jsonl(
    pdb: Path, structure: bts.AtomArray, generated_resids: list[int]
) -> dict[str, dict[str, list[int]]]:
    """
    Given the original unified structure and that with linkers generated,
    determine which residues are in contact with the generated linker and make a .jsonl
    file for ProteinMPNN that would fix all other residues.

    format: {"PDB.stem": {"ChainID": [residue numbers tied position 1]}, {"ChainID: [residue numbers tied position 2]}}
    """
    symmetric_positions_dict = {pdb.stem: []}

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
        symmetric_positions_dict[pdb.stem].append(
            {structure.chain_id[0]: symmetric_positions}
        )

    return symmetric_positions_dict
