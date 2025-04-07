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
    This function detects a synthetic strategy involving late-stage amide formation
    (in the first half of the synthesis) from an ester and an amine.
    """
    amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                ester_found = False
                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(
                            ester_pattern
                        )
                    ):
                        ester_found = True
                        break

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]C")
                amine_found = False
                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(
                            amine_pattern
                        )
                    ):
                        amine_found = True
                        break

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("C(=O)N")
                amide_found = False
                if (
                    product
                    and Chem.MolFromSmiles(product)
                    and Chem.MolFromSmiles(product).HasSubstructMatch(amide_pattern)
                ):
                    amide_found = True

                if ester_found and amine_found and amide_found:
                    print("Late-stage amide formation detected at depth", depth)
                    amide_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return amide_formation_detected
