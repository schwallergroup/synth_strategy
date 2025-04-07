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
    This function detects if the synthesis uses a convergent approach where a
    trifluoromethyl-containing fragment is coupled in a late stage.
    """
    # Track if we've found a late-stage coupling with a CF3-containing fragment
    late_stage_cf3_coupling = False

    # SMARTS pattern for trifluoromethyl group
    cf3_pattern = Chem.MolFromSmarts("[#6]([F])([F])[F]")

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cf3_coupling

        if (
            node["type"] == "reaction"
            and depth <= 1
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if one of the reactants contains a CF3 group
            cf3_reactants = 0
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(cf3_pattern):
                    cf3_reactants += 1
                    print(f"Found CF3-containing reactant at depth {depth}: {reactant}")

            # If we have multiple reactants and at least one contains CF3, it's a convergent coupling
            if len(reactants) >= 2 and cf3_reactants >= 1:
                late_stage_cf3_coupling = True
                print(f"Found late-stage CF3 coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_cf3_coupling
