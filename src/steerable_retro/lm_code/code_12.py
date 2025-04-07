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
    This function detects late-stage C=C bond formation via Heck-type cross-coupling for fragment diversification.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Heck coupling pattern
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and len(reactant_mols) >= 2:
                    # Check for aryl halide in one reactant
                    halide_pattern = Chem.MolFromSmarts("[c]-[I,Br]")
                    # Check for alkene in another reactant
                    alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
                    # Check for styrene-like pattern in product
                    styrene_pattern = Chem.MolFromSmarts("[c]-[C]=[C]")

                    has_halide = any(
                        r and r.HasSubstructMatch(halide_pattern) for r in reactant_mols
                    )
                    has_alkene = any(
                        r and r.HasSubstructMatch(alkene_pattern) for r in reactant_mols
                    )
                    has_styrene = product_mol.HasSubstructMatch(styrene_pattern)

                    if has_halide and has_alkene and has_styrene:
                        print("Found late-stage Heck-type coupling for fragment diversification")
                        found_pattern = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern
