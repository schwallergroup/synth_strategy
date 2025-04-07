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
    Detects a synthetic strategy involving late-stage ether formation through alkylation.
    """
    has_late_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_late_ether_formation

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("metadata", {}).get("depth", -1)

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth > 1:
                return

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if not product_mol or not all(reactant_mols):
                    return

                # Check for ether formation
                ether_pattern = Chem.MolFromSmarts("[#6]O[#6]")
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Br,Cl,I,F]")

                if (
                    product_mol.HasSubstructMatch(ether_pattern)
                    and any(r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols if r)
                    and any(r.HasSubstructMatch(alkyl_halide_pattern) for r in reactant_mols if r)
                ):
                    has_late_ether_formation = True
                    print(f"Detected late-stage ether formation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_late_ether_formation
