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
    Detects if the synthesis route includes an amide formation
    (coupling of carboxylic acid with amine)
    """
    found_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit mol objects
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if product_mol and all(reactant_mols):
                    # SMARTS for carboxylic acid
                    acid_pattern = Chem.MolFromSmarts("[#6][C](=[O])[OH]")

                    # SMARTS for amine
                    amine_pattern = Chem.MolFromSmarts("[NH2][#6]")

                    # SMARTS for amide
                    amide_pattern = Chem.MolFromSmarts("[#6][C](=[O])[NH][#6]")

                    # Check if reactants contain acid and amine
                    has_acid = any(r.HasSubstructMatch(acid_pattern) for r in reactant_mols)
                    has_amine = any(r.HasSubstructMatch(amine_pattern) for r in reactant_mols)

                    # Check if product has amide
                    has_amide = product_mol.HasSubstructMatch(amide_pattern)

                    if has_acid and has_amine and has_amide:
                        print("Found amide formation")
                        found_amide_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_amide_formation
