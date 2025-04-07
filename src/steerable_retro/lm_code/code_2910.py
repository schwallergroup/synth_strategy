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
    This function detects if the synthesis employs a protection-based aldehyde synthesis strategy,
    specifically looking for acetal deprotection to form an aldehyde.
    """
    found_strategy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # In retrosynthetic direction: aldehyde -> acetal
            # In forward direction: acetal -> aldehyde

            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check if any reactant contains an acetal pattern
            acetal_pattern = Chem.MolFromSmarts("[O;!H0]-[C;!H0]-[O;!H0]")
            cyclic_acetal_pattern = Chem.MolFromSmarts("[O;R]-[C;R]-[O;R]")
            aldehyde_pattern = Chem.MolFromSmarts("[C;H1](=O)")

            if product_mol and product_mol.HasSubstructMatch(aldehyde_pattern):
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and (
                        reactant_mol.HasSubstructMatch(acetal_pattern)
                        or reactant_mol.HasSubstructMatch(cyclic_acetal_pattern)
                    ):
                        found_strategy = True
                        print(f"Found acetal deprotection to aldehyde at depth {depth}")
                        break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_strategy
