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
    This function detects a synthetic strategy involving TMS protection and deprotection
    of an alcohol.
    """
    has_tms_protection = False
    has_tms_deprotection = False

    def dfs_traverse(node):
        nonlocal has_tms_protection, has_tms_deprotection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for TMS patterns
                tms_pattern = Chem.MolFromSmarts("[O][Si]([C])([C])[C]")
                alcohol_pattern = Chem.MolFromSmarts("[OH]")

                # Check for TMS deprotection
                reactant_has_tms = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(tms_pattern):
                            reactant_has_tms = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        reactant_has_tms
                        and product_mol
                        and product_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        has_tms_deprotection = True
                        print("Found TMS deprotection")
                except:
                    pass

                # Check for TMS protection
                reactant_has_alcohol = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(alcohol_pattern):
                            reactant_has_alcohol = True
                            break
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        reactant_has_alcohol
                        and product_mol
                        and product_mol.HasSubstructMatch(tms_pattern)
                    ):
                        has_tms_protection = True
                        print("Found TMS protection")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have either protection or deprotection
    strategy_present = has_tms_protection or has_tms_deprotection
    print(f"TMS protection/deprotection strategy detected: {strategy_present}")
    return strategy_present
