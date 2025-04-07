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
    This function detects nitro reduction to amine in the synthetic route.
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

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            nitro_in_reactants = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                    nitro_in_reactants = True
                    break

            # Check for amine group in product where nitro was
            if nitro_in_reactants:
                product_mol = Chem.MolFromSmiles(product)
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                    # This is a simplification - ideally we would check that the NH2 is at the same position
                    # where the nitro group was, but that requires atom mapping
                    nitro_reduction_found = True
                    print("Detected nitro reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_reduction_found
