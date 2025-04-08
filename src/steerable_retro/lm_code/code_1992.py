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
    This function detects if the synthesis route employs a late-stage halogenation strategy,
    specifically looking for halogenation (I, Br, Cl, F) in the final steps.
    """
    late_stage_halogenation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_halogenation_found

        if node["type"] == "reaction" and depth <= 1:  # Check only final two steps
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains halogen that wasn't in reactants
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Check for halogens in product
                halogen_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")
                if product_mol.HasSubstructMatch(halogen_pattern):
                    # Check if this halogen was newly introduced
                    all_reactants_have_halogen = True
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and not reactant_mol.HasSubstructMatch(halogen_pattern):
                            all_reactants_have_halogen = False
                            break

                    if not all_reactants_have_halogen:
                        print(f"Late-stage halogenation detected at depth {depth}")
                        late_stage_halogenation_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_halogenation_found
