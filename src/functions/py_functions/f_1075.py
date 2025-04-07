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
    This function detects if the synthesis route employs ester hydrolysis to form a carboxylic acid.
    """
    # Flag to track if we found ester hydrolysis
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        # Check if this is a reaction node
        if node["type"] == "reaction":
            # Get the reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Create RDKit mol objects
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Check for ester in reactants
                    ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                    # Check for carboxylic acid in product
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                    # Check conditions for ester hydrolysis
                    has_ester = any(
                        mol and mol.HasSubstructMatch(ester_pattern)
                        for mol in reactant_mols
                    )
                    has_carboxylic_acid = product_mol and product_mol.HasSubstructMatch(
                        carboxylic_acid_pattern
                    )

                    if has_ester and has_carboxylic_acid:
                        print("Found ester hydrolysis to form carboxylic acid")
                        found_ester_hydrolysis = True
                except:
                    print("Error processing reaction SMILES")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_ester_hydrolysis
