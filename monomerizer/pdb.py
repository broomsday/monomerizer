"""
Functions for manipulating and analyzing PDB structures.
"""


from pathlib import Path

import numpy as np
import biotite.structure as bts
from biotite.structure.io import load_structure
from volumizer import volumizer
from volumizer import utils as volumizer_utils

from monomerizer.constants import (
    RESNUM_OFFSET,
    BACKBONE_ATOMS,
    MONOMER_CHAIN,
)
from monomerizer.paths import FILL_RES_DIR


def renumber_residues(
    structure: bts.AtomArray, chain_resnum_offset: dict[str, int]
) -> bts.AtomArray:
    """
    For each residue in `structure`, increment the residue number by the value for that chain in `chain_resnum_offset`.
    """
    resids = [atom.res_id + chain_resnum_offset[atom.chain_id] for atom in structure]
    return resids


def unify_chain_ids(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Make all chain-ids the same

    Also renumbers all atoms to start at round multiples of 100.
    """
    # determine the starting resnum for each chain
    chain_ids = []
    for chain_id in structure.chain_id:
        if chain_id not in chain_ids:
            chain_ids.append(chain_id)

    chain_resnum_offset = {
        chain_id: i * RESNUM_OFFSET for i, chain_id in enumerate(chain_ids)
    }

    # rename res_ids and chain_ids
    structure.res_id = renumber_residues(structure, chain_resnum_offset)
    structure.chain_id = [MONOMER_CHAIN] * len(structure)

    return structure


def mutate_to_res(
    structure: bts.AtomArray, fixed_res_ids: list[int], res_name: str = "ALA"
) -> bts.AtomArray:
    """
    Mutate all residues to alanine resname, and delete all sidechains.
    """
    structure = structure[np.isin(structure.atom_name, BACKBONE_ATOMS)]
    structure.res_name = [
        res_name if atom.res_id not in fixed_res_ids else atom.res_name
        for atom in structure
    ]

    return structure


def trim_termini(structure: bts.AtomArray, n_trim: int, c_trim: int) -> bts.AtomArray:
    """
    Remove residues from the termini prior to diffusion.

    This can be useful for helping to bias the directionality of the diffused segments.

    Assumes all chains are composed of the same residues.
    """
    for chain in bts.chain_iter(structure):
        res = chain[np.isin(chain.atom_name, "CA")].res_id
        n_res = res[:n_trim]
        c_res = res[::-1][:c_trim]
        trim_res = np.concatenate([n_res, c_res])
        break

    structure = structure[~np.isin(structure.res_id, trim_res)]

    return structure


def compute_symmetry_rmsd(structure: bts.AtomArray, resids: list[int]) -> float:
    """
    Assuming newly generated residues are symmetric and link previously symmetric residue groups,
    compute the mean RMSD between all symmetry units.
    """
    resid_segments = []
    last_resid = None
    segment = []
    for resid in resids:
        if (last_resid is not None) and (resid - 1 > last_resid):
            resid_segments.append(segment)
            segment = []

        segment.append(resid)
        last_resid = resid
    resid_segments.append(segment)

    structure = structure[structure.atom_name == "CA"]
    structure_segments = [
        structure[np.isin(structure.res_id, resids)] for resids in resid_segments
    ]

    ref_id = int(len(structure_segments) / 2)
    ref_segment = structure_segments[ref_id]

    mobile_segments = [
        segment for id, segment in enumerate(structure_segments) if id != ref_id
    ]
    rmsds = []
    for mobile_segment in mobile_segments:
        superimposed, _ = bts.superimpose(ref_segment, mobile_segment)
        rmsds.append(bts.rmsd(ref_segment, superimposed))

    return np.mean(rmsds)


def get_resids_in_contact(
    structure: bts.AtomArray, query_res_ids: list[int], cutoff_distance: float = 4.0
) -> list[int]:
    """
    Given a structure and query list of residue IDs, return all resids with atoms within a cutoff distance.
    Note: Does not include the query list itself.
    """
    atom_idx_to_resid = {
        atom_idx: atom.res_id for atom_idx, atom in enumerate(structure)
    }
    query_atom_idx = [
        atom_idx
        for atom_idx, res_id in atom_idx_to_resid.items()
        if res_id in query_res_ids
    ]

    cell_list = bts.CellList(structure, cell_size=cutoff_distance)
    adjacency_matrix = cell_list.create_adjacency_matrix(cutoff_distance)

    contacting_res_ids = []
    for atom_idx, res_id in atom_idx_to_resid.items():
        if res_id in query_res_ids:
            continue
        if res_id in contacting_res_ids:
            continue

        if True in adjacency_matrix[atom_idx, query_atom_idx]:
            contacting_res_ids.append(res_id)

    return contacting_res_ids


def get_total_residues(pdb_file: Path) -> int:
    """
    Get the total number of residues in the structure.
    """
    structure = load_structure(pdb_file)
    return len(structure[structure.atom_name == "CA"].res_id)


def fill_pore(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Add alanine residues within the pore of the structure.
    """
    volumizer_utils.set_resolution(5.0)
    _, aligned_structure, volumes_structure = volumizer.volumize_structure(structure)

    # pull out the pore
    pore_voxel_structure = volumes_structure[volumes_structure.res_name == "POR"]

    # setup the resids for downstream ProteinGenerator use
    pore_start_res_id = max(structure.res_id) + RESNUM_OFFSET
    pore_voxel_structure.res_id = [
        pore_start_res_id + i for i in range(len(pore_voxel_structure))
    ]
    pore_res_ids = list(pore_voxel_structure.res_id)

    # load and prepare the residue to be used for pore filling
    fill_res_structure = load_structure(FILL_RES_DIR / "ALA.pdb")
    pore_res_chain_id = chr(ord(aligned_structure.chain_id[0]) + 1)
    fill_res_structure.chain_id = [pore_res_chain_id] * len(fill_res_structure)

    # replace each pore atom with the fill residue
    pore_structure = bts.AtomArray(0)
    for i, pore_voxel in enumerate(pore_voxel_structure):
        displacement = bts.displacement(bts.centroid(fill_res_structure), pore_voxel)
        fill_res_structure = bts.translate(fill_res_structure, displacement)
        fill_res_structure.res_id = [pore_res_ids[i]] * len(fill_res_structure)
        pore_structure += fill_res_structure

    return (pore_res_ids, pore_res_chain_id, (aligned_structure + pore_structure))


def get_buried_res_ids(
    structure: bts.AtomArray, buried_cutoff: float = 10.0
) -> list[int]:
    """
    Return a list of all res-ids that are buried from solvent.

    Considers only the sidechain (and CA) atoms.

    Definition being < 10 A^2 SASA across the sidechain atoms only.
    """
    residue_starts = bts.get_residue_starts(structure)

    sasa_per_atom = bts.sasa(structure)
    sasa_per_residue = [
        np.sum(sasa_per_atom[residue_start : residue_starts[i + 1]])
        if i + 1 < len(residue_starts)
        else np.sum(sasa_per_atom[residue_start:])
        for i, residue_start in enumerate(residue_starts)
    ]

    return [
        res_id
        for i, res_id in enumerate(structure[structure.atom_name == "CA"].res_id)
        if sasa_per_residue[i] <= buried_cutoff
    ]


def clean_generated_structure(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Remove the pore and the C-terminal duplicated atom.
    """
    # remove the pore
    structure = structure[structure.chain_id == MONOMER_CHAIN]

    # remove the last residue which is a duplicate of the first
    last_res_id = structure.res_id[-1]
    structure = structure[structure.res_id != last_res_id]

    return structure
