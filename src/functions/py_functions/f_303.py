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
    This function detects a synthetic strategy involving late-stage deprotection of a tert-butyl ester.
    """
    found_deprotection = False
    depth_of_deprotection = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal found_deprotection, depth_of_deprotection

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for tert-butyl ester deprotection
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                all(mol is not None for mol in reactant_mols)
                and product_mol is not None
            ):
                tbu_ester_pattern = Chem.MolFromSmarts("[C](=O)OC(C)(C)C")
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

                has_tbu_ester = any(
                    len(mol.GetSubstructMatches(tbu_ester_pattern)) > 0
                    for mol in reactant_mols
                )
                has_carboxylic_acid = (
                    len(product_mol.GetSubstructMatches(carboxylic_acid_pattern)) > 0
                )

                if has_tbu_ester and has_carboxylic_acid:
                    found_deprotection = True
                    depth_of_deprotection = min(depth_of_deprotection, depth)
                    print(
                        f"Detected tert-butyl ester deprotection at depth {depth}: {rsmi}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if deprotection is found at a late stage (low depth)
    result = found_deprotection and depth_of_deprotection <= 1
    print(
        f"Late-stage deprotection detected: {result} (depth: {depth_of_deprotection})"
    )
    return result
