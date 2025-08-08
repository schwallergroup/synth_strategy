#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects an early-stage aryl-nitrogen coupling strategy,
    where an aryl group is coupled with a nitrogen heterocycle early in the synthesis.
    """
    found_early_aryl_n_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_early_aryl_n_coupling

        if node["type"] == "reaction" and depth >= 2:  # Early stage (depth >= 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl-nitrogen coupling
                aryl_f_pattern = Chem.MolFromSmarts("[c][F]")
                aryl_nh_pattern = Chem.MolFromSmarts("[c][nH]")
                aryl_n_aryl_pattern = Chem.MolFromSmarts("[c][n][c]")

                has_aryl_f = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(aryl_f_pattern) for r in reactants if r
                )
                has_aryl_nh = any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(aryl_nh_pattern) for r in reactants if r
                )
                has_aryl_n_aryl_product = (
                    Chem.MolFromSmiles(product).HasSubstructMatch(aryl_n_aryl_pattern)
                    if product
                    else False
                )

                if (has_aryl_f or has_aryl_nh) and has_aryl_n_aryl_product:
                    print(f"Found early-stage aryl-nitrogen coupling at depth {depth}")
                    found_early_aryl_n_coupling = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_early_aryl_n_coupling
