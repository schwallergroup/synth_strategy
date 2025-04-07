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
    This function detects a late-stage nitrile hydrolysis to form an amide.
    """
    found_nitrile_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_nitrile_hydrolysis

        if (
            node["type"] == "reaction" and node.get("depth", 0) <= 1
        ):  # Late stage (low depth)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitrile hydrolysis
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactant_mol is not None and product_mol is not None:
                    # Check for nitrile in reactant
                    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                    # Check for amide in product
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                    if (
                        reactant_mol.HasSubstructMatch(nitrile_pattern)
                        and product_mol.HasSubstructMatch(amide_pattern)
                        and not reactant_mol.HasSubstructMatch(amide_pattern)
                    ):
                        print("Found late-stage nitrile hydrolysis")
                        found_nitrile_hydrolysis = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitrile_hydrolysis
