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
    This function detects if the synthetic route involves multiple amide bond formations
    (at least 2 separate reactions forming amide bonds).
    """
    amide_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol:
                    # Amide pattern
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    # Count amide bonds in reactants and product
                    reactant_amide_count = sum(
                        [
                            len(mol.GetSubstructMatches(amide_pattern)) if mol else 0
                            for mol in reactant_mols
                        ]
                    )
                    product_amide_count = (
                        len(product_mol.GetSubstructMatches(amide_pattern)) if product_mol else 0
                    )

                    # If product has more amide bonds than reactants combined, amide formation occurred
                    if product_amide_count > reactant_amide_count:
                        print(f"Found amide formation at depth {depth}")
                        amide_formation_count += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if at least 2 amide formations were detected
    return amide_formation_count >= 2
