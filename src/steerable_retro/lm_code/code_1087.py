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
    This function detects if the synthetic route uses Boc protection/deprotection
    of nitrogen-containing groups.
    """
    has_boc_protection = False

    def dfs_traverse(node):
        nonlocal has_boc_protection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc group in reactants or product
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")

                has_boc_in_reactants = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(boc_pattern):
                            has_boc_in_reactants = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    has_boc_in_product = product_mol and product_mol.HasSubstructMatch(boc_pattern)

                    # If Boc appears or disappears, it's a protection/deprotection
                    if (has_boc_in_reactants and not has_boc_in_product) or (
                        not has_boc_in_reactants and has_boc_in_product
                    ):
                        has_boc_protection = True
                        print(f"Detected Boc protection/deprotection: {rsmi}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_boc_protection
