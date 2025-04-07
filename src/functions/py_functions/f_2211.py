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
    Detects conversion of a nitrile group to a primary amine.
    """
    has_nitrile_to_amine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile_to_amine

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check if reactant has nitrile and product has primary amine
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product)

                        if reactant_mol and product_mol:
                            nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
                            amine_pattern = Chem.MolFromSmarts("[#6]-[#7;H2]")

                            if reactant_mol.HasSubstructMatch(
                                nitrile_pattern
                            ) and product_mol.HasSubstructMatch(amine_pattern):
                                has_nitrile_to_amine = True
                                print(
                                    f"Found nitrile to amine conversion at depth {depth}"
                                )
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_nitrile_to_amine
