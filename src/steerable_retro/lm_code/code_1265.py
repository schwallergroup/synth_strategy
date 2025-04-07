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
    This function detects a strategy involving late-stage esterification of a phenol.
    """
    found_late_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_esterification

        if node["type"] == "reaction" and depth <= 1:  # Late stage = low depth (0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Phenol pattern
                    phenol_pattern = Chem.MolFromSmarts("[OH]-[c]")
                    # Ester pattern
                    ester_pattern = Chem.MolFromSmarts("[c]-[O]-[C](=[O])-[#6]")

                    # Check for phenol in reactants and ester in product
                    if reactant_mol.HasSubstructMatch(
                        phenol_pattern
                    ) and product_mol.HasSubstructMatch(ester_pattern):
                        found_late_esterification = True
                        print(f"Found late-stage phenol esterification at depth {depth}: {rsmi}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage phenol esterification detected: {found_late_esterification}")
    return found_late_esterification
