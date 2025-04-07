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
    Detects if the route contains multiple amide bond formations at different steps.
    """
    amide_formation_depths = []

    def dfs_traverse(node):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node["metadata"].get("depth", -1)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acyl chloride or carboxylic acid in reactants
            acyl_pattern = Chem.MolFromSmarts("[C](=[O])[Cl,OH]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

            has_acyl = False
            has_amine = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(acyl_pattern):
                        has_acyl = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product)

            if (
                has_acyl
                and has_amine
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                print(f"Found amide formation at depth {depth}")
                amide_formation_depths.append(depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 2 amide formations at different depths
    if len(amide_formation_depths) >= 2:
        print(f"Found sequential amide formations at depths: {amide_formation_depths}")
        return True
    return False
