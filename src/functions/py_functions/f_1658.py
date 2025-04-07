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
    Detects if the synthetic route includes a Boc protection/deprotection sequence.
    """
    boc_protected_intermediates = []
    boc_deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_protected_intermediates, boc_deprotection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                boc_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")

                # Check if product has Boc group (protection)
                if product_mol.HasSubstructMatch(boc_pattern):
                    boc_protected_intermediates.append(Chem.MolToSmiles(product_mol))
                    print("Boc protection detected")

                # Check if reactant has Boc group but product doesn't (deprotection)
                if reactant_mol.HasSubstructMatch(
                    boc_pattern
                ) and not product_mol.HasSubstructMatch(boc_pattern):
                    boc_deprotection_found = True
                    print("Boc deprotection detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if both protection and deprotection were found
    return len(boc_protected_intermediates) > 0 and boc_deprotection_found
