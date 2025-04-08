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
    This function detects if the synthesis route employs a protection/deprotection sequence,
    particularly focusing on amine protection with phthalimide.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phthalimide protection (NH2 → phthalimide)
            amine_pattern = Chem.MolFromSmarts("[NH2][#6]")
            phthalimide_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)N1[#6]")

            # Check for protection
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and amine_pattern
                    and phthalimide_pattern
                    and reactant_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(phthalimide_pattern)
                ):
                    print(f"Amine protection detected: {rsmi}")
                    protection_found = True
                    break

            # Check for deprotection
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and amine_pattern
                    and phthalimide_pattern
                    and reactant_mol.HasSubstructMatch(phthalimide_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    print(f"Amine deprotection detected: {rsmi}")
                    deprotection_found = True
                    break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return (
        protection_found or deprotection_found
    )  # Return True if either protection or deprotection is found
