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
    Detects if the synthetic route involves multiple amide bond formations.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                amide_pattern = Chem.MolFromSmarts("[#7][#6](=[#8])")
                matches_in_product = len(product_mol.GetSubstructMatches(amide_pattern))

                # Count amide bonds in reactants
                amide_in_reactants = 0
                for r in reactants:
                    reactant_mol = Chem.MolFromSmiles(r) if r else None
                    if reactant_mol:
                        amide_in_reactants += len(reactant_mol.GetSubstructMatches(amide_pattern))

                # If product has more amide bonds than reactants combined, amide formation occurred
                if matches_in_product > amide_in_reactants:
                    print(
                        f"Found amide formation: {matches_in_product - amide_in_reactants} new amide bonds"
                    )
                    amide_formation_count += matches_in_product - amide_in_reactants

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_count >= 2  # Return True if at least 2 amide formations
