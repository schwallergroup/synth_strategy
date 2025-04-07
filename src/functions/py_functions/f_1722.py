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
    This function detects a synthetic strategy involving a Weinreb amide intermediate
    followed by a late-stage multi-component reaction (MCR).
    """
    # Initialize tracking variables
    has_weinreb_amide = False
    has_mcr = False

    def dfs_traverse(node):
        nonlocal has_weinreb_amide, has_mcr

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for multi-component reaction (4+ components)
            if len(reactants_smiles) >= 4:
                print(f"Found MCR with {len(reactants_smiles)} components")
                has_mcr = True

            # Check for Weinreb amide in reactants
            weinreb_amide_pattern = Chem.MolFromSmarts(
                "[#6]-[#8]-[#7](-[#6])-[#6](=[#8])-[#6]"
            )
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(weinreb_amide_pattern):
                        print("Found Weinreb amide in reactants")
                        has_weinreb_amide = True
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return has_weinreb_amide and has_mcr
