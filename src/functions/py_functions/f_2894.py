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
    This function detects Boc protection of a primary amine in the synthetic sequence.
    """
    boc_protection_used = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_used

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if reactant contains primary amine
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            primary_amine_pattern = Chem.MolFromSmarts("[#7;H2]")
            has_primary_amine = any(
                mol and mol.HasSubstructMatch(primary_amine_pattern)
                for mol in reactant_mols
                if mol
            )

            # Check if product contains Boc-protected amine
            product_mol = Chem.MolFromSmiles(product_smiles)
            boc_pattern = Chem.MolFromSmarts("[#6][#8][#6](=[#8])[#7]")
            has_boc = product_mol and product_mol.HasSubstructMatch(boc_pattern)

            if has_primary_amine and has_boc:
                print(f"Detected Boc protection at depth {depth}")
                boc_protection_used = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return boc_protection_used
