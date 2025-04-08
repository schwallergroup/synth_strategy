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
    This function detects if the synthesis route employs amide formation from a carboxylic acid and amine.
    """
    # Flag to track if we found amide formation
    found_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_amide_formation

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

                    # Check for carboxylic acid in reactants
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                    # Check for amine in reactants
                    amine_pattern = Chem.MolFromSmarts("[NH]")

                    # Check for amide in product
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    # Check conditions for amide formation
                    has_carboxylic_acid = any(
                        mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                        for mol in reactant_mols
                    )
                    has_amine = any(
                        mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                    )
                    has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                    if has_carboxylic_acid and has_amine and has_amide:
                        print("Found amide formation from carboxylic acid and amine")
                        found_amide_formation = True
                except:
                    print("Error processing reaction SMILES")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_amide_formation
