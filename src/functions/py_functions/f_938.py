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
    Detects SNAr reaction for C-N bond formation, typically with
    chloropyridine and amine nucleophile.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for chloropyridine or similar leaving group on heteroaromatic
            halo_heteroaromatic_pattern = Chem.MolFromSmarts("[n]:[c]-[Cl,F,Br]")

            # Check for amine nucleophile
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            # Check for C-N bond formation in product
            cn_bond_pattern = Chem.MolFromSmarts("[n]:[c]-[NH]")

            has_halo_heteroaromatic = False
            has_amine = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(halo_heteroaromatic_pattern):
                        has_halo_heteroaromatic = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product)
            has_cn_bond = False
            if product_mol and product_mol.HasSubstructMatch(cn_bond_pattern):
                has_cn_bond = True

            if has_halo_heteroaromatic and has_amine and has_cn_bond:
                print("SNAr C-N bond formation detected")
                snar_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_detected
