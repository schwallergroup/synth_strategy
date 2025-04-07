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
    This function detects a sequence of nitro reduction to amine followed by isocyanate formation.
    """
    nitro_reduction = False
    amine_to_isocyanate = False
    nitro_reduction_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction, amine_to_isocyanate, nitro_reduction_depth

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[NH2]")
                ):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[N+](=O)[O-]")
                        ):
                            nitro_reduction = True
                            nitro_reduction_depth = depth
                            print(f"Found nitro reduction at depth {depth}")

                # Check for amine to isocyanate
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[O]=[C]=[N]")
                ):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[NH2]")
                        ):
                            amine_to_isocyanate = True
                            # Check if this happens after nitro reduction
                            if (
                                nitro_reduction_depth != -1
                                and depth < nitro_reduction_depth
                            ):
                                print(
                                    f"Found amine to isocyanate at depth {depth}, after nitro reduction"
                                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return true if both transformations are found and in the correct sequence
    return nitro_reduction and amine_to_isocyanate and nitro_reduction_depth != -1
