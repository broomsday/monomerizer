"""
Functions to aid in monomerizing a homo-oligomer.
"""


from pathlib import Path
from dataclasses import dataclass
import subprocess
import random

import numpy as np
from tqdm import tqdm
import biotite.structure as bts
from biotite.structure.io import load_structure, save_structure
import jsonlines

from monomerizer import pdb, proteinmpnn
from monomerizer.constants import (
    ANGSTROMS_PER_HELICAL_RESIDUE,
    PEPTIDE_BOND_ATOMS,
    UPPER_HELICAL_MODIFIER,
    LOWER_HELICAL_MODIFIER,
    EXTENDED_RESIDUES_PER_LOOP,
    LOOP_DISTANCE,
)


random.seed()


@dataclass
class GeneratedStructure:
    """Information on a ProteinGenerator generated structure."""

    pdb_file: Path
    total_length: int
    linker_length: int
    symmetry_rmsd: float

    def __init__(self, input_pdb: Path, generated_pdb: Path, linker_length: int):
        input_structure = load_structure(input_pdb)
        generated_structure = load_structure(generated_pdb)

        # remove the pore and the duplicated final residue
        generated_structure = pdb.clean_generated_structure(generated_structure)
        save_structure(generated_pdb.with_suffix(".clean.pdb"), generated_structure)

        generated_resids = get_unified_resids(
            input_structure, generated_structure, linker_length
        )

        fixed_positions = proteinmpnn.make_fixed_positions_jsonl(
            generated_pdb, generated_structure, generated_resids
        )
        with open(generated_pdb.with_suffix(".fixed.jsonl"), mode="wb") as jsonl_file:
            with jsonlines.Writer(jsonl_file) as writer:
                writer.write(fixed_positions)

        symmetric_positions = proteinmpnn.make_symmetric_positions_jsonl(
            generated_pdb, generated_structure, generated_resids
        )
        with open(
            generated_pdb.with_suffix(".symmetric.jsonl"), mode="wb"
        ) as jsonl_file:
            with jsonlines.Writer(jsonl_file) as writer:
                writer.write(symmetric_positions)

        self.pdb_file = generated_pdb
        self.total_length = pdb.get_total_residues(generated_pdb)
        self.linker_length = linker_length
        self.symmetry_rmsd = pdb.compute_symmetry_rmsd(
            generated_structure, generated_resids
        )


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


def generate_base_contigs(structure: bts.AtomArray) -> str:
    """
    Based on the spacing between N- and C-termini of adjacent chains,
    generates a `contigs` string for use in RFDiffusion that should be
    capable of generating a helical segment to span the termini.
    """
    peptide_bond_atoms = structure[np.isin(structure.atom_name, PEPTIDE_BOND_ATOMS)]

    last_atom = peptide_bond_atoms[0]
    contig_string = f"{peptide_bond_atoms[0].chain_id}{peptide_bond_atoms[0].res_id}-"
    for atom in peptide_bond_atoms:
        # when we hit a chain break, add our notation
        if (atom.res_id - last_atom.res_id) > 1:
            contig_string += f"{last_atom.res_id}"
            contig_string += ",XXX_LINKER,"
            contig_string += f"{atom.chain_id}{atom.res_id}-"
        last_atom = atom

    contig_string += f"{last_atom.res_id}"
    contig_string += ",XXX_LINKER,"
    contig_string += f"{peptide_bond_atoms[0].chain_id}{peptide_bond_atoms[0].res_id}-{peptide_bond_atoms[0].res_id}"

    return contig_string


def parse_monomer_length(contigs: str) -> int:
    """
    Assuming the input file is sequentially res-numbered, determine the length of the monomers.
    """
    lengths = [
        (int(contig.split("-")[-1]) - int(contig.split("-")[0][1:]) + 1)
        for contig in contigs.split(",")
        if "-" in contig
    ]
    lengths.pop()  # remove the terminal duplicate residue

    return int(np.mean(lengths))


def generate_pore_contigs(
    pore_res_ids: list[int],
    pore_res_chain_id: str,
    linker_length: int,
    base_contigs: str,
) -> tuple[str, int]:
    """
    Add contigs for the pore-filling individual amino acids that will be used.

    Additionally and a N->C linker residue.

    Return these contigs, as well as the number of extra symmetry units that have been added.
    """
    monomer_length = parse_monomer_length(base_contigs)
    symmetry_length = monomer_length + linker_length

    symmetry_unit_count = int(len(pore_res_ids) / symmetry_length)
    pore_res_needed = (
        symmetry_unit_count * symmetry_length
    ) - 1  # -1 for duplicated terminal residue

    if len(pore_res_ids) >= pore_res_needed:
        keep_pore_res_ids = random.sample(pore_res_ids, pore_res_needed)
    else:
        keep_pore_res_ids = random.choices(pore_res_ids, pore_res_needed)

    pore_contigs = ",0 "
    for res_id in keep_pore_res_ids:
        pore_contigs += f"{pore_res_chain_id}{res_id}-{res_id},0 "

    return (pore_contigs, symmetry_unit_count)


def make_shell_script(
    input_pdb: Path,
    output_dir: Path,
    res_length: int,
    contigs: str,
    inpaint: str,
    extra_symmetry: int,
    num_designs: int,
    steps: int,
) -> str:
    """
    Write a shell script that will run ProteinGenerator using the input pdb file and computed contigs.
    """
    shell_str = "#!/bin/bash\n\n"
    shell_str += "cd /home/broom/AlphaCarbon/software/protein_generator\n"
    shell_str += "source /usr/etc/profile.d/conda.sh\n"
    shell_str += "conda deactivate\n"
    shell_str += "conda activate proteingenerator\n"

    contig_map = contigs.replace("XXX_LINKER", str(res_length))
    symmetry = contigs.count("XXX_LINKER")

    shell_str += "python ./inference.py "
    shell_str += f"--contigs {contig_map} "
    shell_str += f"--inpaint_seq {inpaint} "
    shell_str += f"--symmetry {symmetry + extra_symmetry} "
    shell_str += f"--num_designs {num_designs} "
    shell_str += f"--out {output_dir.absolute()}/{input_pdb.stem}_{res_length}_monomer "
    shell_str += "--save_best_plddt "
    # shell_str += "--predict_symmetric "
    shell_str += f"--pdb {input_pdb.absolute()} "
    shell_str += f"--T {steps}\n"
    shell_str += "cd /home/broom/AlphaCarbon/code/porepep\n"

    return shell_str


def generate_linked_structures(
    input_pdb: Path,
    pore_res_ids: list[int],
    pore_res_chain_id: str,
    output_dir: Path,
    res_lengths: list[int],
    num_designs: int,
    contigs: str,
    inpaint: str,
    steps: int = 25,
) -> list[GeneratedStructure]:
    """
    Run ProteinGenerator to generate linker residues across a series of lengths.
    """
    scans = []
    for res_length in tqdm(res_lengths, desc="Generating structures"):
        scan_dir = output_dir / str(res_length)
        scan_dir.mkdir(exist_ok=True, parents=True)

        # dont't run if we already have enough generated structures (e.g. from previously stopped run)
        if len(list(scan_dir.glob("*.pdb"))) < num_designs:
            pore_contigs, extra_symmetry = generate_pore_contigs(
                pore_res_ids, pore_res_chain_id, res_length, contigs
            )
            shell_script = make_shell_script(
                input_pdb,
                scan_dir,
                res_length,
                contigs + pore_contigs,
                inpaint,
                extra_symmetry,
                num_designs,
                steps,
            )

            subprocess.run(
                shell_script, shell=True, stdout=subprocess.DEVNULL, check=False
            )

        for design_idx in range(num_designs):
            run_pdb = (
                scan_dir / f"{input_pdb.stem}_{res_length}_monomer_{design_idx:06d}.pdb"
            )
            scans.append(GeneratedStructure(input_pdb, run_pdb, res_length))

    return scans


def get_unified_resids(
    input_structure: bts.AtomArray,
    generated_structure: bts.AtomArray,
    linker_length: int,
) -> list[int]:
    """
    Given the input structure and ProteinGenerator generated structure, return a list of the generated resids.
    """
    # ProteinGenerator renumbers residues to start from 1 and be completely sequential
    #   so adjust resids to match
    input_resids = input_structure[input_structure.atom_name == "CA"].res_id

    unified_resids = []
    last_resid = None
    chain_breaks = 0
    for resid in input_resids:
        if last_resid is None:
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

    # return unified_resids, generated_resids
    return generated_resids


def get_fixed_res_ids(structure: bts.AtomArray, pore_res_ids: list[int]) -> list[int]:
    """
    Determine which residues should be fixed to their input identity:

    1. If they are in contact with the pore
    2. If they have very low solvent exposure
    """
    pore_contacting_res_ids = pdb.get_resids_in_contact(
        structure, pore_res_ids, cutoff_distance=8.0
    )
    buried_res_ids = pdb.get_buried_res_ids(
        structure[~np.isin(structure.res_id, pore_res_ids)]
    )

    return pore_contacting_res_ids + buried_res_ids


def generate_inpaint_contig(
    structure: bts.AtomArray, fixed_res_ids: list[int], pore_chain_id: str
) -> str:
    """
    Given residues
    """
    protein = structure[~np.isin(structure.chain_id, [pore_chain_id])]
    protein_chain_id = protein.chain_id[0]
    res_ids = protein[protein.atom_name == "CA"].res_id
    inpaint_res_ids = sorted(list(set(res_ids).difference(fixed_res_ids)))

    inpaint_contig = ""
    for res_id in inpaint_res_ids:
        inpaint_contig += f"{protein_chain_id}{res_id}-{res_id},"
    inpaint_contig = inpaint_contig.rstrip(",")

    return inpaint_contig
