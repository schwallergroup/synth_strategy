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
    This function detects if the synthesis involves a late-stage esterification
    (formation of an ester bond in the final or penultimate step).
    """
    has_late_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_esterification

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check final or penultimate steps
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if an ester is formed
                reactants_mol = [
                    Chem.MolFromSmiles(r) for r in reactants_part.split(".")
                ]
                product_mol = Chem.MolFromSmiles(product_part)

                # Define ester pattern
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O;!H1]")

                # Check if ester is in product but not in reactants
                if product_mol and product_mol.HasSubstructMatch(ester_pattern):
                    ester_in_reactants = False
                    for r_mol in reactants_mol:
                        if r_mol and r_mol.HasSubstructMatch(ester_pattern):
                            ester_in_reactants = True
                            break

                    if not ester_in_reactants:
                        has_late_esterification = True
                        print(f"Found late-stage esterification at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(
        f"Synthesis {'has' if has_late_esterification else 'does not have'} late-stage esterification"
    )
    return has_late_esterification
