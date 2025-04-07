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
    Detects if the route contains a protection-deprotection sequence,
    specifically looking for carbamate (Cbz) protection and deprotection.
    """
    protection_found = False
    deprotection_found = False
    protection_depth = -1
    deprotection_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal protection_found, deprotection_found, protection_depth, deprotection_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carbamate protection (formation of Cbz group)
                if not protection_found:
                    # Look for carbamate formation: presence in product but not in reactants
                    carbamate_pattern = "[NX3]-[CX3](=[OX1])-[OX2]-[CH2]-[cX3]:[cX3]"
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(carbamate_pattern)
                    ):
                        # Check if reactants don't have this pattern
                        has_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts(carbamate_pattern)
                            ):
                                has_in_reactants = True
                                break

                        if not has_in_reactants:
                            protection_found = True
                            protection_depth = depth
                            print(f"Carbamate protection found at depth {depth}")

                # Check for carbamate deprotection
                if not deprotection_found:
                    # Look for carbamate removal: presence in reactants but not in product
                    carbamate_pattern = "[NX3]-[CX3](=[OX1])-[OX2]-[CH2]-[cX3]:[cX3]"
                    primary_amine_pattern = "[NX3H2]"

                    # Check if product has primary amine
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts(primary_amine_pattern)
                    ):
                        # Check if any reactant has carbamate
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts(carbamate_pattern)
                            ):
                                deprotection_found = True
                                deprotection_depth = depth
                                print(f"Carbamate deprotection found at depth {depth}")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if protection happened before deprotection
    if (
        protection_found
        and deprotection_found
        and protection_depth > deprotection_depth
    ):
        print("Protection-deprotection sequence detected")
        return True

    return False
