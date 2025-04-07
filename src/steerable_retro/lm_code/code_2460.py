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
    Detects benzyl ether protection of phenol in the synthetic route.
    """
    found_benzyl_protection = False

    def dfs_traverse(node):
        nonlocal found_benzyl_protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if not product_mol or not all(reactant_mols):
                return

            # Check for benzyl protection
            phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")
            benzyl_bromide_pattern = Chem.MolFromSmarts("[Br][CH2][c]1[cH][cH][cH][cH][cH]1")
            benzyl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[CH2]-[c]1[cH][cH][cH][cH][cH]1")

            has_phenol = any(mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols)
            has_benzyl_bromide = any(
                mol.HasSubstructMatch(benzyl_bromide_pattern) for mol in reactant_mols
            )
            forms_benzyl_ether = product_mol.HasSubstructMatch(benzyl_ether_pattern)

            if has_phenol and has_benzyl_bromide and forms_benzyl_ether:
                print("Found benzyl ether protection of phenol")
                found_benzyl_protection = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_benzyl_protection
