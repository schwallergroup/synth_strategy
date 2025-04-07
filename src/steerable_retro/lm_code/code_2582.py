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
    This function detects if the synthesis uses a late-stage coupling to connect
    two complex fragments.
    """
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1 and len(reactants) >= 2:
                # Check if both reactants are complex (more than 15 atoms)
                complex_reactants = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 15:
                        complex_reactants += 1

                if complex_reactants >= 2:
                    late_stage_coupling = True
                    print(f"Late-stage fragment coupling detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage fragment coupling strategy detected: {late_stage_coupling}")
    return late_stage_coupling
