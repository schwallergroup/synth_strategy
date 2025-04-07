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
    This function detects if the synthesis includes an ester hydrolysis to carboxylic acid.
    """
    ester_hydrolysis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an ester hydrolysis reaction
                ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
                acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")

                # Check if ester in reactant and acid in product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if (
                        reactant_mol
                        and product_mol
                        and reactant_mol.HasSubstructMatch(ester_pattern)
                        and product_mol.HasSubstructMatch(acid_pattern)
                        and not product_mol.HasSubstructMatch(ester_pattern)
                    ):
                        ester_hydrolysis_found = True
                        print(f"Ester hydrolysis found at depth {depth}")
                        break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return ester_hydrolysis_found
