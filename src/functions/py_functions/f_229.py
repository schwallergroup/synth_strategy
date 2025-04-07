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
    This function detects if the synthetic route involves an early-stage
    bromination followed by a coupling reaction.
    """
    bromination_depth = -1
    coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal bromination_depth, coupling_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for bromination
                    if "Br" in product and any("Br" in r for r in reactants):
                        product_mol = Chem.MolFromSmiles(product)
                        aryl_bromide_pattern = Chem.MolFromSmarts("c-[#35]")

                        if product_mol is not None and product_mol.HasSubstructMatch(
                            aryl_bromide_pattern
                        ):
                            bromination_depth = depth
                            print(f"Bromination detected at depth {depth}")

                    # Check for coupling reaction
                    if len(reactants) >= 2:
                        product_mol = Chem.MolFromSmiles(product)
                        biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

                        if product_mol is not None and product_mol.HasSubstructMatch(
                            biaryl_pattern
                        ):
                            coupling_depth = depth
                            print(f"Coupling reaction detected at depth {depth}")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if bromination occurs at a higher depth (earlier in synthesis) than coupling
    early_bromination = (
        bromination_depth > coupling_depth
        and bromination_depth != -1
        and coupling_depth != -1
    )
    if early_bromination:
        print(
            f"Early-stage bromination (depth {bromination_depth}) followed by coupling (depth {coupling_depth})"
        )

    return early_bromination
