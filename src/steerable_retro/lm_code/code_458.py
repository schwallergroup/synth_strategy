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
    Detects if the synthesis uses Boc protection of amines.
    """
    has_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc anhydride or similar reagent in reactants
                boc_reagent_pattern = re.compile(r"CC\(C\)\(C\)OC\(=O\)O", re.IGNORECASE)

                for reactant in reactants:
                    if boc_reagent_pattern.search(reactant):
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                # Check if product has Boc group
                                boc_group = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                                if product_mol.HasSubstructMatch(boc_group):
                                    has_boc_protection = True
                                    print(f"Found Boc protection at depth {depth}")
                        except:
                            continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_boc_protection
