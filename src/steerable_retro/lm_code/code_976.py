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
    This function detects if the synthesis uses Boc protection of an amine,
    particularly from a nitrile precursor.
    """
    has_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile to Boc-protected amine transformation
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                boc_pattern = Chem.MolFromSmarts("[N]C(=O)OC(C)(C)C")

                has_nitrile = False
                has_boc_reagent = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(nitrile_pattern):
                            has_nitrile = True
                        if "OC(=O)OC(C)(C)C" in reactant:
                            has_boc_reagent = True

                product_mol = Chem.MolFromSmiles(product)
                has_boc_protected = product_mol and product_mol.HasSubstructMatch(boc_pattern)

                if (has_nitrile or has_boc_reagent) and has_boc_protected:
                    print(f"Detected Boc protection strategy at depth {depth}")
                    has_boc_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_boc_protection
