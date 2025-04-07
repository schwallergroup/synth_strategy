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
    Detects if the synthesis involves ester hydrolysis as a key step.
    """
    has_ester_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant has ester and product has carboxylic acid
                has_ester = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                        if mol.HasSubstructMatch(ester_pattern):
                            has_ester = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and has_ester:
                    acid_pattern = Chem.MolFromSmarts("C(=O)O")
                    if product_mol.HasSubstructMatch(acid_pattern):
                        print(f"Found ester hydrolysis at depth {depth}")
                        has_ester_hydrolysis = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_ester_hydrolysis
