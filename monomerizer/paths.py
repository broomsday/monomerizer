"""
Define project paths.
"""

from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT_DIR / "data"
FILL_RES_DIR = DATA_DIR / "fill_residues"

PROTEIN_MPNN_DIR = "/home/broom/AlphaCarbon/software/ProteinMPNN"
