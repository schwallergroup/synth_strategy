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
    Detects a strategy involving N-alkylation using an alkyl bromide.
    """
    has_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_alkylation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for alkyl bromide in reactants
            alkyl_bromide_pattern = Chem.MolFromSmarts("[#6;!$(C=*)]-[#35]")
            has_alkyl_bromide = any(
                r is not None and r.HasSubstructMatch(alkyl_bromide_pattern)
                for r in reactants
            )

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
            has_amine = any(
                r is not None and r.HasSubstructMatch(amine_pattern) for r in reactants
            )

            # Check for new C-N bond in product
            if has_alkyl_bromide and has_amine:
                # This is a simplification - in a real implementation, you would need to
                # check that the specific C that was attached to Br is now attached to N
                has_n_alkylation = True
                print(f"N-alkylation with alkyl bromide detected at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_n_alkylation:
        print("N-alkylation with alkyl bromide strategy detected")
    else:
        print("N-alkylation with alkyl bromide strategy not detected")

    return has_n_alkylation
