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
    This function detects if the route uses thioether (C-S bond) formation
    as a key connection strategy for fragment assembly.
    """
    thioether_formation_detected = False

    def dfs_traverse(node):
        nonlocal thioether_formation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiol pattern in reactants
            thiol_pattern = Chem.MolFromSmarts("[SH]")
            thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

            # Check if reactants contain thiol and product contains thioether
            thiol_in_reactants = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(thiol_pattern):
                        thiol_in_reactants = True
                        break
                except:
                    continue

            # Check if product contains thioether
            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    thiol_in_reactants
                    and product_mol
                    and product_mol.HasSubstructMatch(thioether_pattern)
                ):
                    print("Thioether formation detected")
                    thioether_formation_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thioether_formation_detected
