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
    Detects a synthetic route that includes coupling between heterocycles.
    """
    has_heterocycle_coupling = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Define heterocycle patterns
                pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")
                indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and reactants:
                    # Check if reactants contain different heterocycles
                    heterocycle_types = set()
                    for r in reactants:
                        if r:
                            if r.HasSubstructMatch(pyrimidine_pattern):
                                heterocycle_types.add("pyrimidine")
                            if r.HasSubstructMatch(indole_pattern):
                                heterocycle_types.add("indole")

                    # If we have at least 2 different heterocycle types in reactants
                    if len(heterocycle_types) >= 2:
                        has_heterocycle_coupling = True
                        print(f"Found heterocycle coupling between {', '.join(heterocycle_types)}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle coupling strategy detection result: {has_heterocycle_coupling}")
    return has_heterocycle_coupling
