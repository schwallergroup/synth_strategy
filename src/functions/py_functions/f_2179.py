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
    Detects if the synthesis involves multiple amide bond formations.
    """
    # Count amide formations
    amide_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol:
                        # SMARTS for amide
                        amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")

                        # Count amides in product
                        product_amides = len(
                            product_mol.GetSubstructMatches(amide_pattern)
                        )

                        # Count amides in reactants
                        reactant_amides = 0
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                reactant_amides += len(
                                    reactant_mol.GetSubstructMatches(amide_pattern)
                                )

                        # If product has more amides than reactants, amide formation occurred
                        if product_amides > reactant_amides:
                            amide_formation_count += product_amides - reactant_amides
                            print(f"Found amide formation at depth {depth}")
                except:
                    print(
                        "Error processing SMILES in multiple_amide_formation_strategy"
                    )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return amide_formation_count >= 2  # Return True if at least 2 amide formations
