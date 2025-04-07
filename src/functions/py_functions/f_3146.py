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
    Detects the formation of a Weinreb amide (N-methoxy-N-methylamide) in the late stages of synthesis.
    """
    weinreb_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal weinreb_formation_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check late-stage reactions (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a Weinreb amide
                weinreb_pattern = Chem.MolFromSmarts("[#6][#8][#7][#6](=[#8])")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(weinreb_pattern):
                    # Check if Weinreb amide was introduced in this step
                    weinreb_in_reactants = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(weinreb_pattern):
                                weinreb_in_reactants = True
                                break
                        except:
                            continue

                    if not weinreb_in_reactants:
                        print("Detected late-stage Weinreb amide formation")
                        weinreb_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return weinreb_formation_detected
