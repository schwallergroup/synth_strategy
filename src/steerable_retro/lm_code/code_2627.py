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
    Detects a strategy involving late-stage coupling of two complex fragments.
    """
    has_late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_coupling

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at low depths (late in synthesis)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least two reactants
            if len(reactants_smiles) >= 2:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if reactants are complex (have at least 10 atoms)
                complex_reactants = [
                    r for r in reactants if r is not None and r.GetNumAtoms() >= 10
                ]

                if len(complex_reactants) >= 2:
                    # Check if product has more atoms than any single reactant
                    product_atom_count = product.GetNumAtoms()
                    max_reactant_atom_count = max(
                        r.GetNumAtoms() for r in reactants if r is not None
                    )

                    if product_atom_count > max_reactant_atom_count:
                        has_late_coupling = True
                        print(f"Late-stage fragment coupling detected at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_late_coupling:
        print("Late-stage fragment coupling strategy detected")
    else:
        print("Late-stage fragment coupling strategy not detected")

    return has_late_coupling
