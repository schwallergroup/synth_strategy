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
    Detects if the synthesis includes early oxidation of a methyl group to a carboxylic acid.
    """
    found_methyl_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_methyl_oxidation

        if node["type"] == "reaction" and depth >= 3:  # Early-stage reactions
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for methyl pattern in reactants
                methyl_found = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        methyl_pattern = Chem.MolFromSmarts("[CH3]c")
                        if reactant_mol.HasSubstructMatch(methyl_pattern):
                            methyl_found = True

                # Check for carboxylic acid pattern in product
                product_mol = Chem.MolFromSmiles(product)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")

                if (
                    product_mol
                    and methyl_found
                    and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                ):
                    print(f"Found methyl oxidation at depth {depth}")
                    found_methyl_oxidation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_methyl_oxidation
