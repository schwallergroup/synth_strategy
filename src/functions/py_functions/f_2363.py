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
    This function detects if the synthesis includes a nitro reduction step.
    """
    nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
    amine_pattern = Chem.MolFromSmarts("[#7;H2]")
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if (
                reactant_mol
                and product_mol
                and reactant_mol.HasSubstructMatch(nitro_pattern)
                and product_mol.HasSubstructMatch(amine_pattern)
                and not product_mol.HasSubstructMatch(nitro_pattern)
            ):
                nitro_reduction_found = True
                print("Nitro reduction detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
