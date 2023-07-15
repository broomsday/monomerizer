"""
Functions for manipulating PDBs
"""


import itertools

import numpy as np
import biotite.structure as bts

from monomerizer.constants import (
    RESNUM_OFFSET,
    AXIS_ORDER,
    BACKBONE_ATOMS,
)


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
    Make all chain-ids 'A'

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

    # renumber the residues
    structure.res_id = renumber_residues(structure, chain_resnum_offset)

    # rename all chains to be A
    structure.chain_id = ["A"] * len(structure)

    return structure


def align_to_axis(structure: bts.AtomArray, axis: str) -> bts.AtomArray:
    """
    Align the principal axis to the specified axis.
    """
    return bts.orient_principal_components(structure, order=AXIS_ORDER[axis])


def mutate_to_ala(structure: bts.AtomArray) -> bts.AtomArray:
    """
    Mutate all residues to alanine resname, and delete all sidechains.
    """
    structure.res_name = ["ALA"] * len(structure)
    structure = structure[np.isin(structure.atom_name, BACKBONE_ATOMS)]

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


# TODO: abandon this once we have inpainting working as intended
def revert_to_input(
    wt_structure: bts.AtomArray,
    generated_structure: bts.AtomArray,
    wt_resids: list[int],
    generated_resids: list[int],
) -> bts.AtomArray:
    """
    Given an input structure and generated structure, revert all non-generated residues in the generated structure
    to their input amino acids.

    All sidechain atoms are removed.
    """
    wt_resnames = wt_structure[wt_structure.atom_name == "CA"].res_name
    wt_resid_to_resname = {
        resid: resname for resid, resname in zip(wt_resids, wt_resnames)
    }

    generated_resids = generated_structure[generated_structure.atom_name == "CA"].res_id
    generated_resnames = generated_structure[
        generated_structure.atom_name == "CA"
    ].res_name
    generated_resid_resname = {
        resid: resname
        for resid, resname in zip(generated_resids, generated_resnames)
        if resid not in wt_resids
    }

    resid_to_resname = wt_resid_to_resname | generated_resid_resname
    new_resnames = list(
        itertools.chain.from_iterable(
            [
                [resid_to_resname[resid]] * len(BACKBONE_ATOMS)
                for resid in generated_resids
            ]
        )
    )

    # ignore any sidechains that might have been built in by e.g. ProteinGenerator
    new_structure = generated_structure[
        np.isin(generated_structure.atom_name, BACKBONE_ATOMS)
    ]
    new_structure.res_name = new_resnames

    return new_structure


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
