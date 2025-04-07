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
    This function detects the use of Grignard reagent for C-C bond formation,
    particularly for ketone synthesis.
    """
    grignard_addition_detected = False

    def dfs_traverse(node):
        nonlocal grignard_addition_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Grignard reagent (contains Mg)
                has_grignard = False
                for reactant in reactants:
                    if "[Mg]" in reactant:
                        has_grignard = True
                        break

                # Check for ketone formation
                product_mol = Chem.MolFromSmiles(product)
                ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")

                if (
                    has_grignard
                    and product_mol
                    and product_mol.HasSubstructMatch(ketone_pattern)
                ):
                    grignard_addition_detected = True
                    print("Grignard addition for ketone formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return grignard_addition_detected
