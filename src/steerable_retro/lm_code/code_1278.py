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
    This function detects if the synthesis involves incorporation of a phenyl group.
    """
    phenyl_incorporation_detected = False

    # SMARTS pattern for phenyl group
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")

    def dfs_traverse(node):
        nonlocal phenyl_incorporation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(phenyl_pattern):
                # Check if phenyl is not in all reactants
                phenyl_in_all_reactants = True
                for r_mol in reactants_mols:
                    if r_mol and not r_mol.HasSubstructMatch(phenyl_pattern):
                        phenyl_in_all_reactants = False
                        break

                if not phenyl_in_all_reactants:
                    phenyl_incorporation_detected = True
                    print(f"Phenyl incorporation detected in: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Phenyl incorporation strategy: {phenyl_incorporation_detected}")
    return phenyl_incorporation_detected
