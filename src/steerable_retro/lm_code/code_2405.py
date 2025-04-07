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
    This function detects if a synthetic route employs a strategy where a nitrile group
    is maintained through multiple steps and then converted to an isoxazole in the final step.
    """
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
    isoxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#8][#6][#7]1")

    # Track if we found the pattern
    found_pattern = False

    # Track the depth where nitrile disappears
    nitrile_disappears_at_depth = None

    # Track if isoxazole appears in final product
    isoxazole_in_final_product = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, nitrile_disappears_at_depth, isoxazole_in_final_product

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check if this is the final product (depth 0)
                if depth == 0:
                    if mol.HasSubstructMatch(isoxazole_pattern):
                        isoxazole_in_final_product = True
                        print(f"Found isoxazole in final product: {node['smiles']}")

                # Check for nitrile presence
                has_nitrile = mol.HasSubstructMatch(nitrile_pattern)

                # Process children (reactions)
                for child in node.get("children", []):
                    if child["type"] == "reaction":
                        # Get the product of this reaction (which is the current molecule)
                        product_smiles = node["smiles"]
                        product_mol = mol

                        # Get reactants
                        reactants = []
                        for reactant_node in child.get("children", []):
                            if reactant_node["type"] == "mol":
                                reactant_mol = Chem.MolFromSmiles(reactant_node["smiles"])
                                if reactant_mol:
                                    reactants.append(reactant_mol)

                        # Check if any reactant has nitrile but product doesn't
                        if not has_nitrile:
                            for reactant_mol in reactants:
                                if reactant_mol.HasSubstructMatch(nitrile_pattern):
                                    nitrile_disappears_at_depth = depth
                                    print(f"Nitrile disappears at depth {depth}")
                                    break

                        # Continue traversal
                        dfs_traverse(child, depth + 1)
        else:  # node["type"] == "reaction"
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern: nitrile disappears in the final step (depth 0)
    # and isoxazole appears in the final product
    if nitrile_disappears_at_depth == 0 and isoxazole_in_final_product:
        found_pattern = True
        print("Found late-stage nitrile to isoxazole transformation")

    return found_pattern
