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
    Detects if the synthesis route has no ring formation steps.
    """
    ring_formation_found = False

    def dfs_traverse(node):
        nonlocal ring_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count rings in reactants and product
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and reactant_mols:
                total_reactant_rings = sum(len(Chem.GetSSSR(mol)) for mol in reactant_mols)
                product_rings = len(Chem.GetSSSR(product_mol))

                if product_rings > total_reactant_rings:
                    ring_formation_found = True
                    print("Ring formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return not ring_formation_found
