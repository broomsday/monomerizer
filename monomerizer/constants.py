RESNUM_OFFSET = 200
PEPTIDE_BOND_ATOMS = ["N", "C"]
BACKBONE_ATOMS = ["N", "C", "CA", "O"]
ANGSTROMS_PER_HELICAL_RESIDUE = (
    1.5  # https://en.wikipedia.org/wiki/Alpha_helix#Geometry_and_hydrogen_bonding
)
UPPER_HELICAL_MODIFIER = 1.25
LOWER_HELICAL_MODIFIER = 0.9

ANGSTROMS_PER_EXTENDED_RESIDUE = (
    3.8  # https://en.wikipedia.org/wiki/Beta_sheet#Geometry
)
EXTENDED_RESIDUES_PER_LOOP = 5  # Aron's arbitrary number
LOOP_DISTANCE = ANGSTROMS_PER_EXTENDED_RESIDUE * EXTENDED_RESIDUES_PER_LOOP
MONOMER_CHAIN = "A"
