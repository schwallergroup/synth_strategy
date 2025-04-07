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
    Detects if the synthesis uses a late-stage sulfonylation strategy where a sulfonyl group
    is added to a nitrogen heterocycle in one of the final steps.
    """
    sulfonyl_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
    final_product_has_sulfonyl = False
    sulfonylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_sulfonyl, sulfonylation_depth

        if node["type"] == "mol":
            if depth == 0:  # Final product
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(sulfonyl_pattern):
                    final_product_has_sulfonyl = True

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and reactants_mol:
                product_has_sulfonyl = product_mol.HasSubstructMatch(sulfonyl_pattern)
                reactants_have_sulfonyl = reactants_mol.HasSubstructMatch(
                    sulfonyl_pattern
                )

                # Check if this reaction introduces a sulfonyl group
                if product_has_sulfonyl and not reactants_have_sulfonyl:
                    sulfonylation_depth = depth
                    print(f"Sulfonylation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Strategy is present if final product has sulfonyl group and it was introduced in a late stage (depth ≤ 1)
    strategy_present = (
        final_product_has_sulfonyl
        and sulfonylation_depth is not None
        and sulfonylation_depth <= 1
    )

    print(f"Late-stage sulfonylation strategy detected: {strategy_present}")
    if strategy_present:
        print(f"Sulfonylation occurred at depth {sulfonylation_depth}")

    return strategy_present
