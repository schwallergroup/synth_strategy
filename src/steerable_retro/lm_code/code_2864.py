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
    This function detects a nitration step followed by halogenation in the synthesis.
    """
    nitration_depths = []
    halogenation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal nitration_depths, halogenation_depths

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitration
                nitro_pattern = Chem.MolFromSmarts("[#6][N+](=[O])[O-]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    # Check if nitro group was not in reactants
                    nitro_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                            nitro_in_reactants = True
                            break

                    if not nitro_in_reactants:
                        print(f"Found nitration at depth {depth}")
                        nitration_depths.append(depth)

                # Check for halogenation
                hydroxyl_pattern = Chem.MolFromSmarts("[#6][OH]")
                chloride_pattern = Chem.MolFromSmarts("[#6][Cl]")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(hydroxyl_pattern):
                        if product_mol and product_mol.HasSubstructMatch(chloride_pattern):
                            print(f"Found halogenation at depth {depth}")
                            halogenation_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if nitration is followed by halogenation
    for nitration_depth in nitration_depths:
        for halogenation_depth in halogenation_depths:
            if halogenation_depth < nitration_depth:  # Remember: lower depth = later in synthesis
                print(
                    f"Found nitration at depth {nitration_depth} followed by halogenation at depth {halogenation_depth}"
                )
                return True

    return False
