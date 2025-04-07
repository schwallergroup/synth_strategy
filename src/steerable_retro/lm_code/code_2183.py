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
    Detects a synthetic strategy involving nitro group reduction to amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro to amine reduction
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amine_pattern = Chem.MolFromSmarts("c-[NH2]")
                    if product_mol.HasSubstructMatch(amine_pattern):
                        # Check if any reactant has nitro group
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                nitro_pattern = Chem.MolFromSmarts("c-[N+](=[O])-[O-]")
                                if reactant_mol.HasSubstructMatch(nitro_pattern):
                                    nitro_reduction_detected = True
                                    print("Nitro to amine reduction detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro to amine reduction strategy detected: {nitro_reduction_detected}")
    return nitro_reduction_detected
