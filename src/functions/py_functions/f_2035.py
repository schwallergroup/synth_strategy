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
    Detects if the synthesis includes amide formation between a carboxylic acid and an amine.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid and amine patterns in reactants
                carboxylic_acid_found = False
                amine_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")
                        amine_pattern = Chem.MolFromSmarts("[NH2;!$(NC=O)]")

                        if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            carboxylic_acid_found = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_found = True

                # Check for amide pattern in product
                product_mol = Chem.MolFromSmiles(product)
                amide_pattern = Chem.MolFromSmarts("[C$(C=O)][NH]")

                if (
                    product_mol
                    and product_mol.HasSubstructMatch(amide_pattern)
                    and (carboxylic_acid_found and amine_found)
                ):
                    print(f"Found amide formation at depth {depth}")
                    found_amide_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_amide_formation
