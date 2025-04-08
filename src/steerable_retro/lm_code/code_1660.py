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
    Detects if the synthetic route includes O-demethylation (methoxy to hydroxyl).
    """
    o_demethylation_found = False

    def dfs_traverse(node):
        nonlocal o_demethylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for methoxy in reactant and hydroxyl in product
            reactant_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                methoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6;H3]")
                hydroxyl_pattern = Chem.MolFromSmarts("[#6]-[#8;H]")

                if reactant_mol.HasSubstructMatch(
                    methoxy_pattern
                ) and product_mol.HasSubstructMatch(hydroxyl_pattern):
                    # This is a simplification; a more robust implementation would track atom mappings
                    o_demethylation_found = True
                    print("O-demethylation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return o_demethylation_found
