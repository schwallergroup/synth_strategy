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
    This function detects if the synthesis route includes formation of an amide bond.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid and amine in reactants, amide in product
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]C(=O)[OH]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            amide_pattern = Chem.MolFromSmarts("[#6]C(=O)[NH][#6]")

            try:
                acid_found = False
                amine_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            acid_found = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine_found = True

                if acid_found and amine_found:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                        print("Amide bond formation detected")
                        amide_formation_found = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return amide_formation_found
