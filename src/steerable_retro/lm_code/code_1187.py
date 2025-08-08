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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    max_fragments_per_step = 0

    def dfs_traverse(node):
        nonlocal max_fragments_per_step

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (exclude simple reagents)
            complex_reactants = 0
            for reactant in reactants:
                # Simple reagents typically have fewer atoms
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                    complex_reactants += 1

            max_fragments_per_step = max(max_fragments_per_step, complex_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Linear synthesis typically has at most 2 complex fragments per step
    is_linear = max_fragments_per_step <= 2
    print(f"Maximum complex fragments per step: {max_fragments_per_step}")
    print(f"Synthesis strategy is {'linear' if is_linear else 'convergent'}")

    return is_linear
