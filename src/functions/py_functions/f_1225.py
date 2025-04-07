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
    This function detects a synthetic strategy involving late-stage nitrogen functionalization,
    specifically sulfonylation of indole nitrogen in the final steps.
    """
    late_stage_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_functionalization

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    # Check for indole N-sulfonylation
                    indole_sulfonyl_pattern = Chem.MolFromSmarts(
                        "c1ccc2c(c1)ccn2S(=O)(=O)[#6]"
                    )
                    if product_mol.HasSubstructMatch(indole_sulfonyl_pattern):
                        # Check if sulfonyl group is being added in this step
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        if not any(
                            r and r.HasSubstructMatch(indole_sulfonyl_pattern)
                            for r in reactant_mols
                        ):
                            print(
                                f"Detected late-stage indole N-sulfonylation at depth {depth}"
                            )
                            late_stage_functionalization = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_functionalization
