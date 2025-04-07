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
    Detects if the synthesis route includes a late-stage ester hydrolysis
    (conversion of ester to carboxylic acid in the final steps)
    """
    found_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal found_hydrolysis

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (depth 0-1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for carboxylic acid
                    acid_pattern = Chem.MolFromSmarts("[#6][C](=[O])[OH]")

                    # SMARTS for ester
                    ester_pattern = Chem.MolFromSmarts("[#6][C](=[O])[O][#6]")

                    # Check if product has carboxylic acid
                    has_acid_in_product = product_mol.HasSubstructMatch(acid_pattern)

                    # Check if any reactant has ester
                    has_ester_in_reactants = any(
                        r.HasSubstructMatch(ester_pattern) for r in reactant_mols
                    )

                    if has_acid_in_product and has_ester_in_reactants:
                        print(f"Found late-stage ester hydrolysis at depth {depth}")
                        found_hydrolysis = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_hydrolysis
