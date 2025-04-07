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
    This function detects if the synthesis involves a Grignard addition to a ketone
    to form a tertiary alcohol.
    """
    has_grignard_addition = False

    def dfs_traverse(node, depth=0):
        nonlocal has_grignard_addition

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ketone in reactants
            ketone_found = False
            grignard_found = False

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                # Check for ketone pattern
                if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[#6]")):
                    ketone_found = True

                # Check for Grignard reagent (approximate, looking for Mg-containing species)
                if "Mg" in reactant:
                    grignard_found = True

            # Check for tertiary alcohol in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            tertiary_alcohol_found = False
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[C]([O])([#6])([#6])[#6]")
            ):
                tertiary_alcohol_found = True

            if ketone_found and grignard_found and tertiary_alcohol_found:
                has_grignard_addition = True
                print(f"Grignard addition to ketone detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Grignard addition to ketone: {has_grignard_addition}")
    return has_grignard_addition
