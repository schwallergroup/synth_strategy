#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthetic route follows a convergent synthesis strategy
    by checking for multiple fragment couplings at different stages.
    """
    early_coupling = False
    late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal early_coupling, late_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If there are multiple reactants, it's a coupling reaction
                if len(reactants) > 1:
                    if depth <= 1:  # Late stage (depth 0-1)
                        late_coupling = True
                        print(f"Found late-stage coupling at depth {depth}: {rsmi}")
                    elif depth >= 3:  # Early stage (depth 3+)
                        early_coupling = True
                        print(f"Found early-stage coupling at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Early coupling: {early_coupling}, Late coupling: {late_coupling}")
    return early_coupling and late_coupling
