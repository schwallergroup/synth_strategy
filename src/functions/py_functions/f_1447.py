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
    This function detects a strategy involving halogenation followed by potential
    dehalogenation or use of halogen for coupling reactions.
    """
    halogenation_depths = []
    dehalogenation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_depths, dehalogenation_depths

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit molecules
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                ]

                if product_mol and reactant_mols:
                    # Check for halogenation (focus on iodine)
                    iodo_pattern = Chem.MolFromSmarts("[#6]-[#53]")

                    # Halogenation: product has halogen that reactants don't have
                    if product_mol.HasSubstructMatch(iodo_pattern) and not any(
                        r.HasSubstructMatch(iodo_pattern) for r in reactant_mols if r
                    ):
                        halogenation_depths.append(depth)
                        print(f"Detected halogenation at depth {depth}")

                    # Dehalogenation or halogen utilization: reactants have halogen that product doesn't
                    if any(
                        r.HasSubstructMatch(iodo_pattern) for r in reactant_mols if r
                    ) and not product_mol.HasSubstructMatch(iodo_pattern):
                        dehalogenation_depths.append(depth)
                        print(
                            f"Detected dehalogenation or halogen utilization at depth {depth}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if strategy is present
    strategy_present = len(halogenation_depths) > 0

    print(f"Halogenation depths: {halogenation_depths}")
    print(f"Dehalogenation depths: {dehalogenation_depths}")
    print(f"Strategy detected: {strategy_present}")

    return strategy_present
