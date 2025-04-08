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
    This function detects if the synthetic route involves reduction of a nitro group to an amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains nitro group and product contains amine
            nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
            amine_pattern = Chem.MolFromSmarts("[#7;H2]")

            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                # Check if any reactant has nitro group
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_pattern):
                        print("Nitro reduction detected in reaction")
                        nitro_reduction_found = True
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return nitro_reduction_found
