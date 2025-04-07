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
    Detects the activation of a carboxylic acid to an acid chloride for subsequent coupling.
    """
    has_acid_chloride_formation = False

    def dfs_traverse(node):
        nonlocal has_acid_chloride_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                acid_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")

                    if (
                        reactant_mol
                        and acid_pattern
                        and reactant_mol.HasSubstructMatch(acid_pattern)
                    ):
                        acid_in_reactants = True
                        break

                # Check for acid chloride in product
                if acid_in_reactants:
                    product_mol = Chem.MolFromSmiles(product)
                    acid_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")

                    if (
                        product_mol
                        and acid_chloride_pattern
                        and product_mol.HasSubstructMatch(acid_chloride_pattern)
                    ):
                        has_acid_chloride_formation = True
                        print(
                            f"Detected acid chloride formation at depth {node.get('depth', 'unknown')}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_acid_chloride_formation
