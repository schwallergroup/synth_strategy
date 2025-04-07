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
    This function detects a strategy involving amide formation from an ester and amine.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                # Check for amide formation
                product_mol = Chem.MolFromSmiles(product_part)

                # Look for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Check if reactants include an ester and an amine
                    reactants = reactants_part.split(".")
                    has_ester = False
                    has_amine = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                            amine_pattern = Chem.MolFromSmarts("[NH2]")

                            if reactant_mol.HasSubstructMatch(ester_pattern):
                                has_ester = True
                            if reactant_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_ester and has_amine:
                        found_amide_formation = True
                        print(
                            f"Found amide formation from ester and amine at depth {depth}"
                        )
            except:
                print(f"Error processing reaction at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_amide_formation
