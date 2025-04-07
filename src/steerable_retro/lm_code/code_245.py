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
    This function detects a strategy where a carboxylic acid is converted to an amide
    in the synthetic route.
    """
    has_acid_to_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_acid_to_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid to amide conversion
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Amide pattern
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    if product_mol.HasSubstructMatch(amide_pattern):
                        # Check if any reactant has carboxylic acid
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                                if reactant_mol.HasSubstructMatch(acid_pattern):
                                    print(
                                        f"Found carboxylic acid to amide conversion at depth {depth}"
                                    )
                                    has_acid_to_amide = True
                                    break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_acid_to_amide
