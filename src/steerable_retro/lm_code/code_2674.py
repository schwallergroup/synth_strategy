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
    This function detects the introduction of a methoxy group, particularly
    replacing a chloro group with a methoxy group.
    """
    has_methoxy_introduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_methoxy_introduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for methoxy introduction
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if not product_mol or not all(reactants_mols):
                return

            # Methoxy and chloro SMARTS patterns
            methoxy_pattern = Chem.MolFromSmarts("CO[#6]")
            chloro_pattern = Chem.MolFromSmarts("Cl[#6]")

            # Check if product has methoxy and any reactant has chloro
            methoxy_in_product = product_mol.HasSubstructMatch(methoxy_pattern)
            chloro_in_reactants = any(
                r and r.HasSubstructMatch(chloro_pattern) for r in reactants_mols
            )

            if methoxy_in_product and chloro_in_reactants:
                print(f"Found methoxy introduction at depth {depth}")
                has_methoxy_introduction = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return has_methoxy_introduction
