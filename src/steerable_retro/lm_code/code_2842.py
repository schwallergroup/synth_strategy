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
    This function detects if the synthesis uses a leaving group activation strategy,
    specifically looking for alcohol to mesylate conversion.
    """
    # Initialize tracking variable
    has_alcohol_to_mesylate = False

    def dfs_traverse(node):
        nonlocal has_alcohol_to_mesylate

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for alcohol and mesylate
            alcohol_pattern = Chem.MolFromSmarts("[#8H1]")
            mesylate_pattern = Chem.MolFromSmarts("[#8][S](=[#8])(=[#8])[#6]")

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for alcohol to mesylate conversion
                if (
                    any(m and m.HasSubstructMatch(alcohol_pattern) for m in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(mesylate_pattern)
                ):
                    print("Detected alcohol to mesylate conversion")
                    has_alcohol_to_mesylate = True

            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Leaving group activation strategy detected: {has_alcohol_to_mesylate}")
    return has_alcohol_to_mesylate
