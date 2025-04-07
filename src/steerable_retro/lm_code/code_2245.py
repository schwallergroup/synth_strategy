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
    This function detects if the synthetic route involves late-stage halogenation (depth 0-1).
    """
    late_stage_halogenation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_halogenation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
            and depth <= 1
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains halogen but reactants don't
            halogen_pattern = Chem.MolFromSmarts("c[Cl,Br,I,F]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(halogen_pattern):
                # Check if any reactant has the same halogen pattern
                halogen_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(halogen_pattern):
                        halogen_in_reactants = True
                        break

                if not halogen_in_reactants:
                    print(f"Late-stage halogenation detected at depth {depth}")
                    late_stage_halogenation_found = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_halogenation_found
