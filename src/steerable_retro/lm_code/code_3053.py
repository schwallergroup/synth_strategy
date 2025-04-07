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
    This function detects the use of phthalimide protection/deprotection for primary amine.
    """
    phthalimide_deprotection_found = False

    def dfs_traverse(node):
        nonlocal phthalimide_deprotection_found

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for phthalimide pattern in reactants
            phthalimide_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)N1")
            amine_pattern = Chem.MolFromSmarts("[NX3;H2]")

            phthalimide_found = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(phthalimide_pattern):
                    phthalimide_found = True

            product_mol = Chem.MolFromSmiles(product)
            amine_formed = product_mol and product_mol.HasSubstructMatch(amine_pattern)

            if phthalimide_found and amine_formed:
                phthalimide_deprotection_found = True
                print("Found phthalimide deprotection to reveal primary amine")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return phthalimide_deprotection_found
