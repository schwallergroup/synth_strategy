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
    Detects reductive amination strategy (aldehyde to amine conversion)
    in the synthetic sequence.
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

            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
            amine_pattern = Chem.MolFromSmarts("[#7]([#6])[#6]")  # Secondary or tertiary amine

            # Check if any reactant has aldehyde and product has amine
            aldehyde_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(aldehyde_pattern):
                    aldehyde_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            amine_in_product = product_mol and product_mol.HasSubstructMatch(amine_pattern)

            if aldehyde_in_reactants and amine_in_product:
                found_pattern = True
                print(f"Found reductive amination pattern at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_pattern
