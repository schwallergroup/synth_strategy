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
    Detects if the synthesis involves a nitro group reduction to an amine.
    """
    nitro_to_amine_detected = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check if any reactant contains a nitro group
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if product has an amine group
                    if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                        for reactant_smiles in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                            if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                                print("Nitro to amine reduction detected")
                                nitro_to_amine_detected = True
                                break
                except Exception as e:
                    print(f"Error in nitro to amine detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_to_amine_detected
