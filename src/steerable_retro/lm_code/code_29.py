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
    This function detects if the synthesis involves O-demethylation (methoxy to hydroxyl).
    """
    demethylation_detected = False

    def dfs_traverse(node):
        nonlocal demethylation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains methoxy and product contains hydroxyl
            methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6H3]")
            hydroxyl_pattern = Chem.MolFromSmarts("[#6]-[#8H]")

            product_mol = Chem.MolFromSmiles(product) if product else None

            for reactant in reactants:
                if not reactant:
                    continue
                reactant_mol = Chem.MolFromSmiles(reactant)
                if (
                    reactant_mol
                    and reactant_mol.HasSubstructMatch(methoxy_pattern)
                    and product_mol
                    and product_mol.HasSubstructMatch(hydroxyl_pattern)
                ):
                    print("Detected O-demethylation")
                    demethylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return demethylation_detected
