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
    This function detects if the synthesis involves a late-stage amination
    (replacing O with N in a heterocycle).
    """
    has_late_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amination

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if NH3 or similar nitrogen source is in reactants
            if "NH3" in reactants_part or "[NH2]" in reactants_part or "[NH3]" in reactants_part:
                print(f"Potential amination reagent found at depth {depth}")

                # Check if product has more N atoms than reactants
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product = Chem.MolFromSmiles(product_part)

                if None not in reactants and product is not None:
                    # Count N atoms in heterocycles in reactants
                    reactant_n_count = sum(
                        sum(
                            1
                            for atom in mol.GetAtoms()
                            if atom.GetSymbol() == "N" and atom.IsInRing()
                        )
                        for mol in reactants
                    )

                    # Count N atoms in heterocycles in product
                    product_n_count = sum(
                        1
                        for atom in product.GetAtoms()
                        if atom.GetSymbol() == "N" and atom.IsInRing()
                    )

                    if product_n_count > reactant_n_count:
                        print(f"Late-stage amination detected at depth {depth}")
                        has_late_amination = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_late_amination
