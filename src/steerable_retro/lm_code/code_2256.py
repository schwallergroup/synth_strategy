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
    This function detects the formation of sulfonamide bonds in the synthesis route.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for patterns indicative of sulfonamide formation
                    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[O])(=[O])[Cl]")
                    amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N#*)]")
                    sulfonamide_pattern = Chem.MolFromSmarts(
                        "[N;!$(N=*);!$(N#*)][S](=[O])(=[O])[#6]"
                    )

                    # Check reactants
                    has_sulfonyl_chloride = any(
                        Chem.MolFromSmiles(r) is not None
                        and Chem.MolFromSmiles(r).HasSubstructMatch(sulfonyl_chloride_pattern)
                        for r in reactants
                    )
                    has_amine = any(
                        Chem.MolFromSmiles(r) is not None
                        and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                        for r in reactants
                    )

                    # Check product
                    product_mol = Chem.MolFromSmiles(product)
                    has_sulfonamide = product_mol is not None and product_mol.HasSubstructMatch(
                        sulfonamide_pattern
                    )

                    if (has_sulfonyl_chloride or has_amine) and has_sulfonamide:
                        print(f"Found sulfonamide formation at depth {depth}")
                        found_pattern = True

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_pattern
