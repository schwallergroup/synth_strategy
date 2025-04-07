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
    Detects if the synthesis route maintains a halogen (particularly iodine)
    throughout the synthesis as a potential handle for further functionalization.
    """
    final_product_has_iodine = False
    intermediate_has_iodine = False
    iodine_introduction_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_iodine, intermediate_has_iodine, iodine_introduction_depth

        if node["type"] == "mol" and depth == 0:
            # This is the final product
            mol = Chem.MolFromSmiles(node["smiles"])
            iodine_pattern = Chem.MolFromSmarts("c[IX1]")
            if mol and mol.HasSubstructMatch(iodine_pattern):
                final_product_has_iodine = True
                print("Final product contains iodine")

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            iodine_pattern = Chem.MolFromSmarts("c[IX1]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(iodine_pattern):
                intermediate_has_iodine = True

                # Check if iodine is introduced in this reaction
                reactants_have_iodine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(iodine_pattern):
                        reactants_have_iodine = True
                        break

                if not reactants_have_iodine:
                    iodine_introduction_depth = depth
                    print(f"Iodine introduced at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if iodine is maintained from an intermediate to the final product
    return final_product_has_iodine and intermediate_has_iodine
