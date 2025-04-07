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
    Detects methyl to aldehyde oxidation in the synthetic sequence.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for methyl pattern in reactants and aldehyde in product
            methyl_pattern = Chem.MolFromSmarts("[#6][CH3]")
            aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")

            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol or not product_mol.HasSubstructMatch(aldehyde_pattern):
                # Skip if product doesn't have aldehyde
                pass
            else:
                # Check if any reactant has methyl group
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(methyl_pattern):
                        found_pattern = True
                        print(f"Found methyl to aldehyde oxidation at depth {depth}")
                        break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_pattern
