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
    This function detects a synthetic strategy involving the conversion of an
    acid chloride to a Weinreb amide.
    """
    # Initialize tracking variables
    has_acid_chloride_to_weinreb = False

    def dfs_traverse(node):
        nonlocal has_acid_chloride_to_weinreb

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for acid chloride in reactants
            acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#17]")
            # Check for Weinreb amide in product
            weinreb_amide_pattern = Chem.MolFromSmarts(
                "[#6]-[#8]-[#7](-[#6])-[#6](=[#8])-[#6]"
            )

            has_acid_chloride = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acid_chloride_pattern):
                        has_acid_chloride = True
                        break
                except:
                    continue

            has_weinreb = False
            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(weinreb_amide_pattern):
                    has_weinreb = True
            except:
                pass

            if has_acid_chloride and has_weinreb:
                print("Found acid chloride to Weinreb amide transformation")
                has_acid_chloride_to_weinreb = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_acid_chloride_to_weinreb
