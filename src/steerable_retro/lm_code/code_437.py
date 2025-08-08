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
    This function detects if the synthesis follows a convergent approach where
    two or more complex fragments are joined in the final steps (depth 0-1).
    """
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple complex reactants (>20 atoms) at depth 0-1,
            # it's likely a convergent synthesis
            if depth <= 1 and len(reactants) >= 2:
                complex_reactants = 0
                for r in reactants:
                    if r.count("c") + r.count("C") > 10:  # Simple heuristic for complexity
                        complex_reactants += 1

                if complex_reactants >= 2:
                    is_convergent = True
                    print(
                        f"Found convergent synthesis at depth {depth} with {complex_reactants} complex fragments"
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_convergent
