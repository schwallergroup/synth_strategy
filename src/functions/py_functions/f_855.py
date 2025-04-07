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
    This function detects a synthetic strategy involving activation of an alcohol
    to a chloride (OH → Cl) as a leaving group.
    """
    alcohol_to_chloride_detected = False

    def dfs_traverse(node):
        nonlocal alcohol_to_chloride_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to chloride conversion
            alcohol_pattern = Chem.MolFromSmarts("[#8H1][#6]")
            chloride_pattern = Chem.MolFromSmarts("[Cl][#6]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(
                    r and r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols
                )
                and product_mol.HasSubstructMatch(chloride_pattern)
            ):
                # Additional check to confirm it's an alcohol to chloride conversion
                # Look for OH disappearance and Cl appearance
                if any(
                    r and r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols
                ) and not any(
                    r and r.HasSubstructMatch(chloride_pattern) for r in reactant_mols
                ):
                    print(f"Detected alcohol to chloride conversion: {rsmi}")
                    alcohol_to_chloride_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Alcohol activation to chloride strategy detected: {alcohol_to_chloride_detected}"
    )
    return alcohol_to_chloride_detected
