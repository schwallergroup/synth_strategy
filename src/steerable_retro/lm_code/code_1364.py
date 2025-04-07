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
    This function detects if the synthesis involves multiple C-N bond formations
    (e.g., amide and urea).
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])[#6]")
                    urea_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])[#7]")

                    # Count matches in product
                    product_amide_matches = len(product_mol.GetSubstructMatches(amide_pattern))
                    product_urea_matches = len(product_mol.GetSubstructMatches(urea_pattern))

                    # Count matches in reactants
                    reactant_amide_matches = 0
                    reactant_urea_matches = 0

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            reactant_amide_matches += len(
                                reactant_mol.GetSubstructMatches(amide_pattern)
                            )
                            reactant_urea_matches += len(
                                reactant_mol.GetSubstructMatches(urea_pattern)
                            )

                    # If product has more amide or urea groups than reactants combined, a C-N bond was formed
                    if (
                        product_amide_matches > reactant_amide_matches
                        or product_urea_matches > reactant_urea_matches
                    ):
                        cn_bond_formations += 1
                        print(f"Detected C-N bond formation, total count: {cn_bond_formations}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    return cn_bond_formations >= 2  # Return True if at least 2 C-N bond formations
