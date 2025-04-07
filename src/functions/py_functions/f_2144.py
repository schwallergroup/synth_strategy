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
    Detects if the synthetic route employs Boc protection of an amine
    """
    boc_protection_detected = False

    def dfs_traverse(node):
        nonlocal boc_protection_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for free amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            # Check for Boc-protected amine in product
            boc_pattern = Chem.MolFromSmarts(
                "[#7]-[C](=[O])-[O]-[C]([CH3])([CH3])[CH3]"
            )

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if (
                product_mol
                and any(
                    mol and mol.HasSubstructMatch(amine_pattern)
                    for mol in reactant_mols
                )
                and product_mol.HasSubstructMatch(boc_pattern)
            ):
                print("Boc protection detected")
                boc_protection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection_detected
