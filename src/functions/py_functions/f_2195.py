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
    Detects if the synthesis uses multiple instances of O-alkylation
    (hydroxyl to methoxy conversion).
    """
    o_alkylation_count = 0

    def dfs_traverse(node):
        nonlocal o_alkylation_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for hydroxyl to methoxy conversion
            hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
            methoxy_pattern = Chem.MolFromSmarts("[O][CH3]")

            has_hydroxyl_reactant = False
            has_methoxy_product = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(hydroxyl_pattern):
                        has_hydroxyl_reactant = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(methoxy_pattern):
                    has_methoxy_product = True
            except:
                pass

            if has_hydroxyl_reactant and has_methoxy_product:
                print("Found O-alkylation (hydroxyl to methoxy)")
                o_alkylation_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Consider it a strategy if there are at least 2 instances
    multiple_instances = o_alkylation_count >= 2

    if multiple_instances:
        print(
            f"Detected multiple O-alkylation strategy ({o_alkylation_count} instances)"
        )

    return multiple_instances
