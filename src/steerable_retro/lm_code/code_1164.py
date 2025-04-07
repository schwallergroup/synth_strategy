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
    This function detects if a nucleophilic aromatic substitution (SNAr) reaction
    is present in the synthesis route, typically involving a fluoronitrobenzene.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for fluoronitrobenzene pattern in reactants
            fluoro_nitro_pattern = Chem.MolFromSmarts("Fc1c([N+](=O)[O-])cccc1")

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(fluoro_nitro_pattern):
                    # Check if the product has a new C-N bond where the F was
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # This is a simplified check - a more robust implementation would
                        # track the specific atoms involved in the substitution
                        if not product_mol.HasSubstructMatch(fluoro_nitro_pattern):
                            snar_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"SNAr reaction detected: {snar_detected}")

    return snar_detected
