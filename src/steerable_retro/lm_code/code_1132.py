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
    This function detects a late-stage fragment coupling strategy where two complex fragments
    are combined in the final steps of the synthesis.
    """
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (low depth)
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Check if we're combining two complex fragments
            if len(reactants) >= 2:
                complex_fragments = 0
                for reactant in reactants:
                    # Define complexity: contains at least one ring or has more than 10 atoms
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            ring_info = mol.GetRingInfo()
                            if ring_info.NumRings() > 0 or mol.GetNumAtoms() > 10:
                                complex_fragments += 1
                    except:
                        pass

                if complex_fragments >= 2:
                    late_stage_coupling = True
                    print(f"Found late-stage coupling at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_coupling
