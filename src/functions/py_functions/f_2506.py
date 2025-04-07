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
    Detects if the synthetic route contains multiple nitro reduction steps.
    """
    nitro_reduction_count = 0

    def dfs_traverse(node):
        nonlocal nitro_reduction_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a nitro reduction reaction
            reactant_mol = Chem.MolFromSmiles(".".join(reactants))
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                reactant_nitro_matches = len(
                    reactant_mol.GetSubstructMatches(nitro_pattern)
                )
                product_nitro_matches = len(
                    product_mol.GetSubstructMatches(nitro_pattern)
                )

                reactant_amine_matches = len(
                    reactant_mol.GetSubstructMatches(amine_pattern)
                )
                product_amine_matches = len(
                    product_mol.GetSubstructMatches(amine_pattern)
                )

                # If nitro groups decreased and amine groups increased, it's likely a reduction
                if (
                    reactant_nitro_matches > product_nitro_matches
                    and product_amine_matches > reactant_amine_matches
                ):
                    nitro_reduction_count += 1
                    print(
                        f"Detected nitro reduction at depth {node.get('depth', 'unknown')}"
                    )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if multiple nitro reductions were detected
    result = nitro_reduction_count >= 2
    print(
        f"Multiple nitro reductions strategy detected: {result} (count: {nitro_reduction_count})"
    )
    return result
