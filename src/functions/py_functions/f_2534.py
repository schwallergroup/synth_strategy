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
    Detects if the synthesis route involves demethylation of a methoxy group to a hydroxyl group.
    """
    demethylation_detected = False

    def dfs_traverse(node):
        nonlocal demethylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for methoxy in reactants
                    methoxy_pattern = Chem.MolFromSmarts("[#8][#6]")
                    hydroxyl_pattern = Chem.MolFromSmarts("[#8H]")

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(hydroxyl_pattern):
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if (
                                mol
                                and mol.HasSubstructMatch(methoxy_pattern)
                                and not mol.HasSubstructMatch(hydroxyl_pattern)
                            ):
                                demethylation_detected = True
                                print(
                                    f"Demethylation detected: {reactant} -> {product}"
                                )
                                break
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Demethylation strategy detected: {demethylation_detected}")
    return demethylation_detected
