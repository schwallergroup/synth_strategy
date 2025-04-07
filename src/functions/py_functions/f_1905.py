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
    This function detects a synthetic strategy involving multiple aryl ether formations
    (at least 3) in a linear synthesis.
    """
    aryl_ether_formation_count = 0

    def dfs_traverse(node):
        nonlocal aryl_ether_formation_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an aryl ether formation
            # Look for phenol in reactants and aryl ether in product
            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
            aryl_ether_pattern = Chem.MolFromSmarts("[c][O][C]")

            reactant_has_phenol = False
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(phenol_pattern):
                        reactant_has_phenol = True
                        break
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                product_has_aryl_ether = product_mol and product_mol.HasSubstructMatch(
                    aryl_ether_pattern
                )
            except:
                product_has_aryl_ether = False

            if reactant_has_phenol and product_has_aryl_ether:
                aryl_ether_formation_count += 1
                print(f"Found aryl ether formation: {rsmi}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Total aryl ether formations: {aryl_ether_formation_count}")
    return aryl_ether_formation_count >= 3
