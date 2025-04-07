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
    This function detects if the synthesis route includes a phthalimide deprotection step
    to reveal a primary amine.
    """
    phthalimide_deprotection_found = False

    def dfs_traverse(node):
        nonlocal phthalimide_deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phthalimide pattern in reactant
            phthalimide_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)N1[#6]")
            amine_pattern = Chem.MolFromSmarts("[NH2][#6]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and phthalimide_pattern
                    and amine_pattern
                    and reactant_mol.HasSubstructMatch(phthalimide_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    print(f"Phthalimide deprotection detected: {rsmi}")
                    phthalimide_deprotection_found = True
                    break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return phthalimide_deprotection_found
