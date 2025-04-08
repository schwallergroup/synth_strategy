#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the final product contains multiple (3+) aromatic/heteroaromatic rings.
    """
    # Get the root node (final product)
    if route["type"] != "mol":
        print("Root node is not a molecule")
        return False

    product_smiles = route["smiles"]
    print(f"Analyzing product: {product_smiles}")

    product_mol = Chem.MolFromSmiles(product_smiles)

    if not product_mol:
        print("Could not parse product molecule")
        return False

    # Count aromatic rings
    aromatic_rings = 0
    ring_info = product_mol.GetRingInfo()

    # Get all rings as atom indices
    atom_rings = ring_info.AtomRings()

    for ring in atom_rings:
        # Check if ring is aromatic (all atoms in the ring are aromatic)
        is_aromatic = all(product_mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if is_aromatic:
            aromatic_rings += 1
            print(f"Found aromatic ring: {[idx for idx in ring]}")

    print(f"Found {aromatic_rings} aromatic rings in final product")
    return aromatic_rings >= 3
