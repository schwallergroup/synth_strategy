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
    Detects a strategy that includes an alkene reduction step.
    """
    has_alkene_reduction = False

    def dfs_traverse(node):
        nonlocal has_alkene_reduction

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # For alkene reduction, we expect one reactant
            if len(reactants_smiles) == 1:
                reactant = Chem.MolFromSmiles(reactants_smiles[0])
                product = Chem.MolFromSmiles(product_smiles)

                if reactant and product:
                    alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")
                    reactant_alkenes = len(reactant.GetSubstructMatches(alkene_pattern))
                    product_alkenes = len(product.GetSubstructMatches(alkene_pattern))

                    # If reactant has more alkenes than product, it's a reduction
                    if reactant_alkenes > product_alkenes:
                        has_alkene_reduction = True
                        print("Alkene reduction detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Has alkene reduction strategy: {has_alkene_reduction}")
    return has_alkene_reduction
