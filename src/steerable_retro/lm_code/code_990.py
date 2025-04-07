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
    This function detects if the synthetic route contains aromatic nitration.
    """
    found_nitration = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitration

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in product but not in reactants
                nitro_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    # Check if nitro group is not present in reactants
                    nitro_in_reactants = False
                    for reactant in reactants:
                        if "N(=O)(O)" in reactant or "[N+](=[O])[O-]" in reactant:
                            continue  # Skip nitrating agent

                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                            nitro_in_reactants = True
                            break

                    if not nitro_in_reactants:
                        print(f"Found aromatic nitration at depth {depth}")
                        found_nitration = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_nitration
