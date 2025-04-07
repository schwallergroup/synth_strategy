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
    Detects if the synthesis uses a late-stage Suzuki coupling to install a heterocycle.
    Specifically looks for C-C bond formation between an aryl bromide and a boronic acid derivative
    in the final or penultimate step.
    """
    late_stage_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_suzuki

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if one reactant contains Br attached to aromatic carbon
                has_aryl_bromide = False
                has_boronic_acid = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for aryl bromide
                        aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
                        if mol.HasSubstructMatch(aryl_bromide_pattern):
                            has_aryl_bromide = True

                        # Check for boronic acid or derivative
                        boronic_pattern = Chem.MolFromSmarts("[#6]B([OX2])[OX2]")
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic_acid = True

                # Check if product has new C-C bond where Br was
                if has_aryl_bromide and has_boronic_acid:
                    print(f"Found potential Suzuki coupling at depth {depth}")
                    late_stage_suzuki = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_suzuki
