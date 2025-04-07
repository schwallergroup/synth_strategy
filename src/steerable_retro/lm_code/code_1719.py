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
    Detects if the synthesis route involves a halogen to nitrile conversion
    on an aromatic or heterocyclic system.
    """
    found_conversion = False

    def dfs_traverse(node):
        nonlocal found_conversion

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Halogen on aromatic/heterocyclic pattern
                halo_pattern = Chem.MolFromSmarts("[c,n]-[Cl,Br,I,F]")
                # Nitrile pattern
                nitrile_pattern = Chem.MolFromSmarts("[c,n]-[C]#[N]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".")]
                product = Chem.MolFromSmiles(product_str)

                # Check if reactants contain halogenated aromatic/heterocycle
                has_halo = any(
                    r is not None and r.HasSubstructMatch(halo_pattern) for r in reactants
                )

                # Check if product contains nitrile on aromatic/heterocycle
                has_nitrile = product is not None and product.HasSubstructMatch(nitrile_pattern)

                if has_halo and has_nitrile:
                    print("Found halogen to nitrile conversion")
                    found_conversion = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_conversion
