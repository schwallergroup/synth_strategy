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
    Detects protection/deprotection sequences (benzyl, Boc) in the synthesis.
    """
    found_benzyl_protection = False
    found_boc_protection = False

    def dfs_traverse(node):
        nonlocal found_benzyl_protection, found_boc_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for benzyl protection
            phenol_pattern = Chem.MolFromSmarts("[OH][c]")
            benzyl_bromide_pattern = Chem.MolFromSmarts("Br[CH2][c]1[cH][cH][cH][cH][cH]1")
            benzyl_ether_pattern = Chem.MolFromSmarts("[O][CH2][c]1[cH][cH][cH][cH][cH]1")

            # Patterns for Boc protection/deprotection
            boc_group_pattern = Chem.MolFromSmarts("[C](=O)O[C]([CH3])([CH3])[CH3]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and all(reactant_mols):
                # Check for benzyl protection
                has_phenol = any(m and m.HasSubstructMatch(phenol_pattern) for m in reactant_mols)
                has_benzyl_bromide = any(
                    m and m.HasSubstructMatch(benzyl_bromide_pattern) for m in reactant_mols
                )
                forms_benzyl_ether = product_mol.HasSubstructMatch(benzyl_ether_pattern)

                if has_phenol and (has_benzyl_bromide or forms_benzyl_ether):
                    print("Found benzyl protection of phenol")
                    found_benzyl_protection = True

                # Check for Boc protection/deprotection
                has_boc = any(m and m.HasSubstructMatch(boc_group_pattern) for m in reactant_mols)
                product_has_boc = product_mol.HasSubstructMatch(boc_group_pattern)

                if has_boc or product_has_boc:
                    print("Found Boc protection/deprotection")
                    found_boc_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_benzyl_protection and found_boc_protection
