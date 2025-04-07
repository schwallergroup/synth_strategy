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
    This function detects if the synthesis route includes a nitro reduction step (NO2 → NH2).
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactant
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                product_mol = Chem.MolFromSmiles(product)

                if (
                    reactant_mol
                    and product_mol
                    and nitro_pattern
                    and amine_pattern
                    and reactant_mol.HasSubstructMatch(nitro_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                ):
                    print(f"Nitro reduction detected: {rsmi}")
                    nitro_reduction_found = True
                    break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return nitro_reduction_found
