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
    This function detects a synthetic strategy involving an early Claisen condensation
    to form a 1,3-diketone.
    """
    early_claisen_found = False

    def dfs_traverse(node, depth=0):
        nonlocal early_claisen_found

        if node["type"] == "reaction" and depth >= 3:  # Early stage = high depth in retrosynthesis
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product is a 1,3-diketone
                product_mol = Chem.MolFromSmiles(product)
                diketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#6]-[#6](=[O])-[#6]")

                if (
                    product_mol
                    and diketone_pattern
                    and product_mol.HasSubstructMatch(diketone_pattern)
                ):
                    # Check if reactants contain ester and ketone
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#8]-[#6]")
                    ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#6]")

                    has_ester = False
                    has_ketone = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if ester_pattern and reactant_mol.HasSubstructMatch(ester_pattern):
                                has_ester = True
                            if ketone_pattern and reactant_mol.HasSubstructMatch(ketone_pattern):
                                has_ketone = True

                    if has_ester and has_ketone:
                        print("Found early Claisen condensation to form 1,3-diketone")
                        early_claisen_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return early_claisen_found
