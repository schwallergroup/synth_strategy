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
    This function detects if the route contains a Michael addition to an α,β-unsaturated nitrile.
    """
    has_michael_addition = False

    def dfs_traverse(node):
        nonlocal has_michael_addition

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants) and len(reactants) >= 2:
                # Check for Michael addition to α,β-unsaturated nitrile
                unsaturated_nitrile_pattern = Chem.MolFromSmarts("[C]=[C][C]#[N]")
                secondary_amine_pattern = Chem.MolFromSmarts("[N;H1]([C])[C]")
                tertiary_amine_pattern = Chem.MolFromSmarts("[N;H0]([C])([C])[C]")

                reactants_have_unsaturated_nitrile = any(
                    r.HasSubstructMatch(unsaturated_nitrile_pattern) for r in reactants if r
                )
                reactants_have_secondary_amine = any(
                    r.HasSubstructMatch(secondary_amine_pattern) for r in reactants if r
                )
                product_has_tertiary_amine = product.HasSubstructMatch(tertiary_amine_pattern)

                if (
                    reactants_have_unsaturated_nitrile
                    and reactants_have_secondary_amine
                    and product_has_tertiary_amine
                ):
                    has_michael_addition = True
                    print(
                        f"Michael addition to α,β-unsaturated nitrile detected in reaction: {rsmi}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Michael addition to α,β-unsaturated nitrile: {has_michael_addition}")
    return has_michael_addition
