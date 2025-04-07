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
    This function detects if the route contains both amide formation and amide reduction steps.
    """
    has_amide_formation = False
    has_amide_reduction = False

    def dfs_traverse(node):
        nonlocal has_amide_formation, has_amide_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for amide formation
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                amine_pattern = Chem.MolFromSmarts("[N]")
                carboxyl_pattern = Chem.MolFromSmarts("[C](=[O])[O,Cl,Br,I,F]")

                product_has_amide = product.HasSubstructMatch(amide_pattern)
                reactants_have_amide = any(
                    r.HasSubstructMatch(amide_pattern) for r in reactants if r
                )
                reactants_have_amine = any(
                    r.HasSubstructMatch(amine_pattern) for r in reactants if r
                )
                reactants_have_carboxyl = any(
                    r.HasSubstructMatch(carboxyl_pattern) for r in reactants if r
                )

                # Amide formation: carboxyl + amine → amide
                if (
                    product_has_amide
                    and not reactants_have_amide
                    and reactants_have_amine
                    and reactants_have_carboxyl
                ):
                    has_amide_formation = True
                    print(f"Amide formation detected in reaction: {rsmi}")

                # Amide reduction: amide → amine
                amine_ch2n_pattern = Chem.MolFromSmarts("[C][N]")
                if (
                    reactants_have_amide
                    and not product_has_amide
                    and product.HasSubstructMatch(amine_ch2n_pattern)
                ):
                    has_amide_reduction = True
                    print(f"Amide reduction detected in reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Bidirectional amide transformations: {has_amide_formation and has_amide_reduction}")
    return has_amide_formation and has_amide_reduction
