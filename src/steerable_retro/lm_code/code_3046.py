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
    This function detects if the synthetic route employs an ester hydrolysis step
    to generate a carboxylic acid for subsequent reactions.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant includes a methyl/ethyl ester
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                # Check if product has a carboxylic acid
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                ester_found = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(ester_pattern):
                            ester_found = True
                            print("Found ester reactant:", reactant)
                    except:
                        continue

                # Check product for carboxylic acid
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        print("Found carboxylic acid in product:", product)
                        if ester_found:
                            result = True
                            print("Detected ester hydrolysis strategy")
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return result
