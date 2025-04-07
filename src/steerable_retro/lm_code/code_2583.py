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
    This function detects if the synthesis uses nucleophilic aromatic substitution
    with a fluoroaromatic compound.
    """
    snAr_found = False

    def dfs_traverse(node):
        nonlocal snAr_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for fluoroaromatic compound
            fluoro_aromatic_pattern = Chem.MolFromSmarts("[F][c]")
            amine_pattern = Chem.MolFromSmarts("[NH]")

            # Check if any reactant is a fluoroaromatic compound
            fluoro_aromatic_present = False
            amine_present = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(fluoro_aromatic_pattern):
                        fluoro_aromatic_present = True
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_present = True

            # Check if product has C-N bond where F was
            product_mol = Chem.MolFromSmiles(product)

            if (
                fluoro_aromatic_present
                and amine_present
                and product_mol
                and not product_mol.HasSubstructMatch(fluoro_aromatic_pattern)
            ):
                snAr_found = True
                print("Nucleophilic aromatic substitution detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Nucleophilic aromatic substitution strategy detected: {snAr_found}")
    return snAr_found
