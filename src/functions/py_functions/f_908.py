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
    Detects if the synthesis involves strategic protection of functional groups.
    """
    has_protection_step = False

    def dfs_traverse(node, depth=0):
        nonlocal has_protection_step

        if node["type"] == "reaction":
            # Check if this is a reaction node
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection patterns
                # 1. Amine to amide protection
                amine_pattern = Chem.MolFromSmarts("[NH2]-[c]")
                amide_pattern = Chem.MolFromSmarts("[NH]-[C](=[O])-[C]")

                # 2. Alcohol protection (TMS)
                alcohol_pattern = Chem.MolFromSmarts("[OH]-[C]")
                protected_alcohol_pattern = Chem.MolFromSmarts("[O]-[C]-[C]-[Si]")

                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for amine protection
                    has_amine_in_reactants = False
                    for reactant in reactants:
                        try:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                amine_pattern
                            ):
                                has_amine_in_reactants = True
                        except:
                            continue

                    if (
                        has_amine_in_reactants
                        and product_mol
                        and product_mol.HasSubstructMatch(amide_pattern)
                    ):
                        has_protection_step = True
                        print(f"Detected amine protection at depth {depth}")

                    # Check for alcohol protection
                    has_alcohol_in_reactants = False
                    for reactant in reactants:
                        try:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                alcohol_pattern
                            ):
                                has_alcohol_in_reactants = True
                        except:
                            continue

                    if (
                        has_alcohol_in_reactants
                        and product_mol
                        and product_mol.HasSubstructMatch(protected_alcohol_pattern)
                    ):
                        has_protection_step = True
                        print(f"Detected alcohol protection at depth {depth}")

                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_protection_step
