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
    Detects a synthetic strategy involving installation of a tosylate or mesylate leaving group.
    """
    leaving_group_installation_detected = False

    def dfs_traverse(node):
        nonlocal leaving_group_installation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tosylation/mesylation pattern
            alcohol_pattern = Chem.MolFromSmarts("[OH]C")
            tosylate_pattern = Chem.MolFromSmarts("OS(=O)(=O)C")

            # Check if reactants contain alcohol and product contains tosylate
            reactant_has_alcohol = False
            for reactant in reactants:
                if Chem.MolFromSmiles(reactant) and Chem.MolFromSmiles(reactant).HasSubstructMatch(
                    alcohol_pattern
                ):
                    reactant_has_alcohol = True
                    break

            product_has_tosylate = False
            if (
                product
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(tosylate_pattern)
            ):
                product_has_tosylate = True

            if reactant_has_alcohol and product_has_tosylate:
                leaving_group_installation_detected = True
                print("Leaving group installation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return leaving_group_installation_detected
