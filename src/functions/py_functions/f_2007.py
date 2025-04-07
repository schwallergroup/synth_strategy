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
    Detects if the synthetic route contains a carboxylic acid protection-deprotection sequence.
    Specifically looks for COOH → COOMe → COOH pattern.
    """
    # Track if we've seen protection and deprotection
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for esterification (protection)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH3]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        carboxylic_acid_pattern
                    ):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            methyl_ester_pattern
                        ):
                            protection_found = True
                            print("Found carboxylic acid protection")

                # Check for ester hydrolysis (deprotection)
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        methyl_ester_pattern
                    ):
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol and product_mol.HasSubstructMatch(
                            carboxylic_acid_pattern
                        ):
                            deprotection_found = True
                            print("Found carboxylic acid deprotection")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
