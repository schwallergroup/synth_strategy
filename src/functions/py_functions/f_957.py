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
    Detects if the synthetic route employs reductive amination for C-N bond formation.
    """
    reductive_amination_found = False

    # SMARTS patterns
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for aldehyde in reactants
                aldehyde_present = False
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.HasSubstructMatch(aldehyde_pattern):
                        aldehyde_present = True
                        break

                # Check for amine in reactants
                amine_present = False
                for r in reactants_smiles:
                    if "[NH2]" in r or "NH2" in r:
                        amine_present = True
                        break

                # If both aldehyde and amine are present, and product doesn't have aldehyde
                # but has a new C-N bond, it's likely reductive amination
                if aldehyde_present and amine_present:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and not product_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ):
                        reductive_amination_found = True
                        print("Detected reductive amination")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return reductive_amination_found
