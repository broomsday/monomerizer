"""
Define project paths.
"""

from pathlib import Path
import os


ROOT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT_DIR / "data"
FILL_RES_DIR = DATA_DIR / "fill_residues"


def get_protein_generator_path() -> Path | None:
    """
    Use an environment variable to determine where ProteinGenerator is installed.
    """
    try:
        return os.environ["PROTEIN_GENERATOR"]
    except KeyError:
        print("You need to install 'protein_generator'")
        print(
            "Then set the 'PROTEIN_GENERATOR' environment variable to point to it's location"
        )
        quit()


def get_conda_source_path() -> Path | None:
    """
    Use an environment variable to determine where the conda source shell script is located.
    """
    try:
        return os.environ["CONDA_SOURCE"]
    except KeyError:
        print(
            "You need to set the 'CONDA_SOURCE' environment variable to point to conda's sourch.sh script"
        )
        quit()
