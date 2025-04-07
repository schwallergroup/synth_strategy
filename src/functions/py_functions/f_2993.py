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
    Detects if the synthesis uses a Wittig reaction for carbon chain extension.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Wittig reaction patterns
            phosphorus_pattern = Chem.MolFromSmarts("[P]")
            aldehyde_pattern = Chem.MolFromSmarts("[#6H](=O)")
            alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check if product contains new alkene
                if product_mol and product_mol.HasSubstructMatch(alkene_pattern):
                    # Check if reactants contain phosphorus and aldehyde
                    has_phosphorus = False
                    has_aldehyde = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(phosphorus_pattern):
                                has_phosphorus = True
                            if reactant_mol.HasSubstructMatch(aldehyde_pattern):
                                has_aldehyde = True

                    if has_phosphorus and has_aldehyde:
                        result = True
                        print("Detected Wittig olefination for carbon chain extension")
            except:
                print("Error processing SMILES in wittig_olefination_strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
