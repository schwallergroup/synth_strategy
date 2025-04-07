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
    This function detects a linear synthesis with carbonyl functional group progression
    (ester → acid → amide → aldehyde → oxime)
    """
    # Define SMARTS patterns for carbonyl functional groups
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
    weinreb_amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[N][OX2][#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")

    # Track which functional groups are found and at what depth
    found_groups = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(ester_pattern):
                    found_groups["ester"] = depth
                if mol.HasSubstructMatch(acid_pattern):
                    found_groups["acid"] = depth
                if mol.HasSubstructMatch(weinreb_amide_pattern):
                    found_groups["weinreb"] = depth
                if mol.HasSubstructMatch(aldehyde_pattern):
                    found_groups["aldehyde"] = depth
                if mol.HasSubstructMatch(oxime_pattern):
                    found_groups["oxime"] = depth

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 3 of the expected functional groups in the correct order
    # Remember: higher depth = earlier in synthesis (retrosynthetic direction)
    required_groups = ["ester", "acid", "weinreb", "aldehyde", "oxime"]
    found_required = [group for group in required_groups if group in found_groups]

    if len(found_required) >= 3:
        # Check if the order is correct (higher depth should come first in synthesis)
        for i in range(len(found_required) - 1):
            current = found_required[i]
            next_group = found_required[i + 1]
            if found_groups[current] <= found_groups[next_group]:
                print(
                    f"Carbonyl progression order violation: {current} at depth {found_groups[current]} comes before {next_group} at depth {found_groups[next_group]}"
                )
                return False

        print(f"Found carbonyl progression strategy with groups: {found_required}")
        return True

    print("Did not find sufficient evidence of carbonyl progression strategy")
    return False
