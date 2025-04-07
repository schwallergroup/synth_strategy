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
    This function detects a synthetic strategy involving sequential protection of nitrogen atoms,
    specifically Boc protection followed by sulfonylation.
    """
    boc_protection_step = False
    sulfonylation_step = False
    boc_depth = -1
    sulfonyl_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_step, sulfonylation_step, boc_depth, sulfonyl_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    # Check for Boc protection
                    boc_pattern = Chem.MolFromSmarts("[#7]-C(=O)OC(C)(C)C")
                    if product_mol.HasSubstructMatch(boc_pattern):
                        # Check if Boc group is being added in this step
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        if not any(r and r.HasSubstructMatch(boc_pattern) for r in reactant_mols):
                            print(f"Detected Boc protection at depth {depth}")
                            boc_protection_step = True
                            boc_depth = depth

                    # Check for sulfonylation
                    sulfonyl_pattern = Chem.MolFromSmarts("[#7]-S(=O)(=O)-[#6]")
                    if product_mol.HasSubstructMatch(sulfonyl_pattern):
                        # Check if sulfonyl group is being added in this step
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        if not any(
                            r and r.HasSubstructMatch(sulfonyl_pattern) for r in reactant_mols
                        ):
                            print(f"Detected sulfonylation at depth {depth}")
                            sulfonylation_step = True
                            sulfonyl_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both protections occurred and in the correct order (sulfonylation after Boc)
    sequential_protection = (
        boc_protection_step and sulfonylation_step and sulfonyl_depth < boc_depth
    )
    if sequential_protection:
        print("Detected sequential nitrogen protection strategy: Boc followed by sulfonylation")

    return sequential_protection
