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
    Detects if the synthesis route involves a reduction pathway from oxime to amine.
    """
    oxime_reduction_found = False
    oxime_depths = []
    amine_formation_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal oxime_reduction_found, oxime_depths, amine_formation_depths

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for oxime pattern in reactants
            oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")

            # Check for amine formation
            amine_pattern = Chem.MolFromSmarts("[CX4][NX3]")

            product_mol = Chem.MolFromSmiles(product)

            # Check if reactants contain oxime and product contains new amine
            reactants_have_oxime = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(oxime_pattern):
                    reactants_have_oxime = True
                    oxime_depths.append(depth)
                    print(f"Oxime detected in reactant at depth {depth}")
                    break

            if (
                reactants_have_oxime
                and product_mol
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                amine_formation_depths.append(depth)
                print(f"Amine formation from oxime detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have both oxime and subsequent amine formation
    if oxime_depths and amine_formation_depths:
        # Check if there's a sequential relationship
        for oxime_depth in oxime_depths:
            for amine_depth in amine_formation_depths:
                if amine_depth <= oxime_depth:  # Remember lower depth means later in synthesis
                    oxime_reduction_found = True
                    print(
                        f"Oxime to amine reduction pathway detected: oxime at depth {oxime_depth}, amine at depth {amine_depth}"
                    )

    return oxime_reduction_found
