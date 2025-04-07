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
    This function detects if the synthetic route involves an early-stage nitration
    of an aromatic ring.
    """
    # Track if we found nitration in early stages
    found_early_nitration = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_early_nitration, max_depth

        # Track maximum depth to determine what's "early stage"
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitration (addition of nitro group)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
                    if product_mol.HasSubstructMatch(nitro_pattern):
                        # Check if nitro group is new (not in reactants)
                        nitro_count_product = len(
                            product_mol.GetSubstructMatches(nitro_pattern)
                        )

                        total_nitro_count_reactants = 0
                        for reactant in reactants:
                            if "N(=O)(O)" in reactant or "[N+](=[O])[O-]" in reactant:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    total_nitro_count_reactants += len(
                                        reactant_mol.GetSubstructMatches(nitro_pattern)
                                    )

                        if nitro_count_product > total_nitro_count_reactants:
                            print(f"Found nitration at depth {depth}")
                            # We'll determine if it's "early stage" after traversal
                            if depth >= 2:  # Assuming depth >= 2 is "early stage"
                                found_early_nitration = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Determine if nitration occurred in early stage (upper half of synthesis)
    return found_early_nitration
