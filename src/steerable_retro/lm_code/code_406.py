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
    Detects a strategy involving a sequence of functional group interconversions,
    such as alkyne reduction and azide reduction.
    """
    found_alkyne_reduction = False
    found_azide_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alkyne_reduction, found_azide_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check for alkyne reduction
                alkyne_pattern = Chem.MolFromSmarts("[C]#[C]")

                try:
                    reactant_mol = Chem.MolFromSmiles(reactants_str)
                    product_mol = Chem.MolFromSmiles(product_str)

                    if (
                        reactant_mol
                        and product_mol
                        and reactant_mol.HasSubstructMatch(alkyne_pattern)
                        and not product_mol.HasSubstructMatch(alkyne_pattern)
                    ):
                        found_alkyne_reduction = True
                        print(f"Found alkyne reduction at depth {depth}")
                except:
                    pass

                # Check for azide reduction
                azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=[N]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                try:
                    for reactant in reactants_str.split("."):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(azide_pattern):
                            product_mol = Chem.MolFromSmiles(product_str)
                            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                                found_azide_reduction = True
                                print(f"Found azide reduction at depth {depth}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if we found multiple functional group interconversions
    return found_alkyne_reduction and found_azide_reduction
