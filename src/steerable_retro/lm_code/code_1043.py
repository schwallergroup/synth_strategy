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
    This function detects if the synthetic route involves early-stage bromination
    (typically at depth >= 3).
    """
    bromination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal bromination_detected

        if node["type"] == "reaction" and depth >= 3:  # Early stage reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains bromine but reactants don't all have it
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and "Br" in product_smiles:
                    # Check if any reactant doesn't have bromine
                    reactant_list = reactants_smiles.split(".")
                    all_reactants_have_bromine = True
                    for r in reactant_list:
                        if "Br" not in r:
                            all_reactants_have_bromine = False
                            break

                    if not all_reactants_have_bromine:
                        print(f"Early-stage bromination detected at depth {depth}")
                        bromination_detected = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return bromination_detected
