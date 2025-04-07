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
    This function detects if the synthetic route involves Boc protection/deprotection sequences.
    """
    boc_protection_count = 0

    def dfs_traverse(node):
        nonlocal boc_protection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for Boc protection pattern
            if Chem.MolFromSmiles(reactants) and Chem.MolFromSmiles(product):
                boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=O)[#7]")
                amine_pattern = Chem.MolFromSmarts("[NH]")

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                # Boc protection: amine -> Boc-protected amine
                if (
                    reactant_mol.HasSubstructMatch(amine_pattern)
                    and not reactant_mol.HasSubstructMatch(boc_pattern)
                    and product_mol.HasSubstructMatch(boc_pattern)
                ):
                    print("Detected Boc protection")
                    boc_protection_count += 1

                # Boc deprotection: Boc-protected amine -> amine
                if (
                    reactant_mol.HasSubstructMatch(boc_pattern)
                    and not product_mol.HasSubstructMatch(boc_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    print("Detected Boc deprotection")
                    boc_protection_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return boc_protection_count >= 2
