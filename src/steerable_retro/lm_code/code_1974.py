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
    This function detects a synthetic strategy involving the formation of an imidazolidinedione
    ring system from simpler precursors.
    """
    # Flag to track if we found imidazolidinedione formation
    found_formation = False

    def is_imidazolidinedione_formation(reactants_smiles, product_smiles):
        """Check if a reaction forms an imidazolidinedione ring"""
        try:
            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if not product or not all(reactants):
                return False

            # Check if product contains imidazolidinedione
            imidazolidinedione_pattern = Chem.MolFromSmarts(
                "[#7]-1-[#6](=[#8])-[#7]-[#6](=[#8])-[#6]-1"
            )
            if len(product.GetSubstructMatches(imidazolidinedione_pattern)) == 0:
                return False

            # Check if none of the reactants contain imidazolidinedione
            for r in reactants:
                if len(r.GetSubstructMatches(imidazolidinedione_pattern)) > 0:
                    return False

            # If we get here, the reaction forms an imidazolidinedione
            return True

        except Exception as e:
            print(f"Error in is_imidazolidinedione_formation: {e}")
            return False

    def dfs_traverse(node):
        nonlocal found_formation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this reaction forms an imidazolidinedione
                if is_imidazolidinedione_formation(reactants, product):
                    found_formation = True
                    print(f"Found imidazolidinedione formation: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Imidazolidinedione formation strategy detected: {found_formation}")
    return found_formation
