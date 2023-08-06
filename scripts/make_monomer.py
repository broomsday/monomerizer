"""
Modifies a PDB file to be used for motif-scaffolding in order to make it monomeric.

1. Changes all chain-ids to be A (so that the protein is a monomer)
2. Aligns the axis of symmetry along the z-axis (requirement of RFDiffusion)
"""


from pathlib import Path
from shutil import copy as copy_file

import typer
import pandas as pd
from biotite.structure.io import load_structure, save_structure

from monomerizer import monomerizer, pdb


def main(
    input_pdb: Path = typer.Argument(..., help="PDB to be modified"),
    output_dir: Path = typer.Argument(..., help="Directory to hold outputs"),
    generator_steps: int = typer.Option(
        default=25, help="How many denoising steps to take."
    ),
    num_designs: int = typer.Option(
        default=3,
        help="How many structures RFDiffusion should generate per linker length",
    ),
    linker_lengths: str = typer.Option(
        default=None,
        help="List of lengths to use for linking, e.g. '10-20' or '12,15,17`, if not provided then choose automatically",
    ),
    trim_n: int = typer.Option(
        default=0, help="How many N-terminal residues to trim from each chain"
    ),
    trim_c: int = typer.Option(
        default=0, help="How many C-terminal residues to trim from each chain"
    ),
    num_selected: int = typer.Option(
        default=3, help="Top designs by symmetry to copy into the top output directory."
    ),
):
    """
    Change all chain-ids to A and align axis of cyclic symmetry to the z-axis for a given input PDB.

    Examples:
    python scripts/make_monomer.py tests/pdbs/2XQS_short.pdb tmp/ --generator-steps=25 --linker-lengths=15 --num-designs=5
    python scripts/make_monomer.py tests/pdbs/2XQS_short.pdb tmp/ --generator-steps=10 --linker-lengths=10-20 --num-designs=3
    python scripts/make_monomer.py tests/pdbs/4JPP_short.pdb tmp/ --generator-steps=15 --linker-lengths=30 --num-designs=1
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    # load, modify, and save the pdb structure
    unified_structure = load_structure(input_pdb)

    # perform structure modifications
    if (trim_n > 0) or (trim_c > 0):
        unified_structure = pdb.trim_termini(unified_structure, trim_n, trim_c)
    unified_structure = pdb.unify_chain_ids(unified_structure)

    # generate a template for what the contig map will look like
    base_contigs = monomerizer.generate_base_contigs(unified_structure)

    # add the pore
    pore_res_ids, pore_res_chain_id, unified_structure = pdb.fill_pore(
        unified_structure
    )
    # determine which residues the pore contacts or are buried, these should be fixed
    fixed_res_ids = monomerizer.get_fixed_res_ids(unified_structure, pore_res_ids)
    inpaint_contigs = monomerizer.generate_inpaint_contig(
        unified_structure, fixed_res_ids, pore_res_chain_id
    )

    unified_pdb = (output_dir / input_pdb.stem).with_suffix(".unified.pdb")
    save_structure(unified_pdb, unified_structure)

    # get the number of lengths/runs we need to complete
    if linker_lengths is None:
        res_lengths = monomerizer.get_res_range(unified_structure)
    else:
        res_lengths = monomerizer.parse_manual_res_range(linker_lengths)

    # scan over residue lengths and generate structures
    structures = monomerizer.generate_linked_structures(
        unified_pdb,
        pore_res_ids,
        pore_res_chain_id,
        output_dir,
        res_lengths,
        num_designs,
        base_contigs,
        inpaint_contigs,
        steps=generator_steps,
    )

    results = {
        "pdb": [structure.pdb_file for structure in structures],
        "linker_length": [structure.linker_length for structure in structures],
        "symmetry_rmsd": [structure.symmetry_rmsd for structure in structures],
    }
    results_df = pd.DataFrame().from_dict(results)
    results_df = results_df.sort_values(by="symmetry_rmsd")
    selected_dir = output_dir / "top"
    selected_dir.mkdir(exist_ok=True)
    for idx in range(min(num_selected, num_designs * len(res_lengths))):
        copy_file(
            results_df.iloc[idx].pdb,
            selected_dir / results_df.iloc[idx].pdb.name,
        )
    results_df.to_csv(output_dir / "results.csv")
    print(results_df)


if __name__ == "__main__":
    typer.run(main)
