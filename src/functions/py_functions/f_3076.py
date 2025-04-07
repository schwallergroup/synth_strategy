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
    Detects a synthetic route with late-stage coupling of two complex fragments.
    """
    has_late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late-stage (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]

                # Count number of significant fragments (complex molecules)
                reactants = reactants_str.split(".")
                significant_fragments = 0

                for r in reactants:
                    if r:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            # Consider a fragment significant if it has at least 10 atoms
                            if mol.GetNumAtoms() >= 10:
                                significant_fragments += 1

                if significant_fragments >= 2:
                    has_late_stage_coupling = True
                    print(
                        f"Found late-stage coupling of {significant_fragments} complex fragments at depth {depth}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage fragment coupling strategy detection result: {has_late_stage_coupling}"
    )
    return has_late_stage_coupling
