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
    This function detects a strategy involving amide bond formation/cleavage
    in a molecule containing an azide functional group.
    """
    has_azide = False
    has_amide_bond_operation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_azide, has_amide_bond_operation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for azide group
                azide_pattern = "[#7-]=[#7+]=[#7]"
                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts(azide_pattern))
                    for r in reactants
                    if r
                ):
                    has_azide = True
                    print(f"Found azide group at depth {depth}")

                # Check for amide bond formation/cleavage
                amide_pattern = "C(=O)[NH]"
                product_mol = Chem.MolFromSmiles(product) if product else None
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for amide in product but not in reactants (formation)
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts(amide_pattern)):
                    if not any(
                        m and m.HasSubstructMatch(Chem.MolFromSmarts(amide_pattern))
                        for m in reactant_mols
                        if m
                    ):
                        has_amide_bond_operation = True
                        print(f"Found amide bond formation at depth {depth}")

                # Check for amide in reactants but not in product (cleavage)
                if any(
                    m and m.HasSubstructMatch(Chem.MolFromSmarts(amide_pattern))
                    for m in reactant_mols
                    if m
                ):
                    if not (
                        product_mol
                        and product_mol.HasSubstructMatch(Chem.MolFromSmarts(amide_pattern))
                    ):
                        has_amide_bond_operation = True
                        print(f"Found amide bond cleavage at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    strategy_present = has_azide and has_amide_bond_operation
    print(f"Azide-containing amide formation/cleavage strategy detected: {strategy_present}")
    return strategy_present
