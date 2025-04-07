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
    This function detects the use of Weinreb amide for controlled reduction to aldehyde
    """
    # Define SMARTS patterns
    weinreb_amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[N][OX2][#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")

    # Track if we find the patterns and at what depth
    found_weinreb_depth = None
    found_aldehyde_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_weinreb_depth, found_aldehyde_depth

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(weinreb_amide_pattern):
                    found_weinreb_depth = depth
                if mol.HasSubstructMatch(aldehyde_pattern):
                    found_aldehyde_depth = depth

        # Check for Weinreb amide to aldehyde transformation in reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants[0]) if reactants else None
            product_mol = Chem.MolFromSmiles(product) if product else None

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(
                    weinreb_amide_pattern
                ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                    found_weinreb_depth = depth
                    found_aldehyde_depth = depth - 1  # Aldehyde appears one step later in synthesis

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both Weinreb amide and aldehyde in the correct sequence
    if found_weinreb_depth is not None and found_aldehyde_depth is not None:
        if found_weinreb_depth > found_aldehyde_depth:  # Weinreb comes before aldehyde in synthesis
            print(
                f"Found Weinreb amide strategy: Weinreb at depth {found_weinreb_depth}, aldehyde at depth {found_aldehyde_depth}"
            )
            return True

    print("Did not find evidence of Weinreb amide strategy")
    return False
