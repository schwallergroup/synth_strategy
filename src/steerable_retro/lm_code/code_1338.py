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
    Detects the use of acid chloride formation from carboxylic acid
    as part of the overall synthetic strategy.
    """
    found_acid_chloride_formation = False

    def dfs_traverse(node):
        nonlocal found_acid_chloride_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride formation
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if (
                    product_mol
                    and any(
                        m and m.HasSubstructMatch(carboxylic_acid_pattern) for m in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(acid_chloride_pattern)
                ):
                    print("Found acid chloride formation")
                    found_acid_chloride_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_acid_chloride_formation
