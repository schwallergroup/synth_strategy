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
    This function detects a protection-deprotection sequence using Boc group.
    """
    boc_protection = False
    boc_deprotection = False

    def dfs_traverse(node):
        nonlocal boc_protection, boc_deprotection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection (primary amine -> Boc-protected amine)
                boc_group_pattern = Chem.MolFromSmarts("[#6]C([#8])([#8][#6]([#6])([#6])[#6])")
                primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
                protected_amine_pattern = Chem.MolFromSmarts("[NH]C(=O)O[C]([CH3])([CH3])[CH3]")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(protected_amine_pattern):
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                primary_amine_pattern
                            ):
                                boc_protection = True
                                print(
                                    "Boc protection detected: Primary amine -> Boc-protected amine"
                                )
                except:
                    pass

                # Check for Boc deprotection (Boc-protected amine -> primary amine)
                try:
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(protected_amine_pattern):
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and product_mol.HasSubstructMatch(primary_amine_pattern):
                                boc_deprotection = True
                                print(
                                    "Boc deprotection detected: Boc-protected amine -> primary amine"
                                )
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection and boc_deprotection
