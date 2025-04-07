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
    This function detects a strategy involving reduction of unsaturation (C=C to C-C)
    during the synthesis.
    """
    reduction_detected = False

    def dfs_traverse(node):
        nonlocal reduction_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_smiles if r
                    ]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(reactant_mols) and product_mol:
                        # Check for C=C pattern in reactants
                        alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")

                        # Count alkenes in reactants and product
                        reactant_alkene_count = sum(
                            len(mol.GetSubstructMatches(alkene_pattern))
                            for mol in reactant_mols
                        )
                        product_alkene_count = len(
                            product_mol.GetSubstructMatches(alkene_pattern)
                        )

                        # If reactants have more alkenes than product, a reduction likely occurred
                        if reactant_alkene_count > product_alkene_count:
                            print(f"Detected reduction of C=C to C-C")
                            reduction_detected = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return reduction_detected
