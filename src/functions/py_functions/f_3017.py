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
    Detects if the synthesis involves formation of a diaryl ether (Ar-O-Ar) linkage.
    """
    found_diaryl_ether = False

    def dfs_traverse(node):
        nonlocal found_diaryl_ether

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Try to detect diaryl ether formation
            try:
                product_mol = Chem.MolFromSmiles(product_str)

                # Check if product contains diaryl ether pattern
                diaryl_ether_pattern = Chem.MolFromSmarts("c[OX2]c")

                if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    # Check reactants to see if they don't have the pattern
                    reactants = reactants_str.split(".")
                    has_pattern_in_reactants = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            diaryl_ether_pattern
                        ):
                            has_pattern_in_reactants = True
                            break

                    # If pattern is in product but not in reactants, it was formed in this step
                    if not has_pattern_in_reactants:
                        print("Found diaryl ether formation")
                        found_diaryl_ether = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_diaryl_ether
