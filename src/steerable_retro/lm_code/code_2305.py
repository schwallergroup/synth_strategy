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
    This function detects a linear synthesis strategy involving protected aniline derivatives.
    """
    steps_with_protected_aniline = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_protected_aniline, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]

                reactant_mol = Chem.MolFromSmiles(reactants)

                if reactant_mol:
                    # Protected aniline pattern (carbamate)
                    protected_aniline_pattern = Chem.MolFromSmarts("[c]-[N]-[C](=[O])-[O]")

                    if reactant_mol.HasSubstructMatch(protected_aniline_pattern):
                        steps_with_protected_aniline += 1
                        print(f"Found protected aniline in step: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a linear synthesis with protected anilines
    is_linear_with_protected_aniline = steps_with_protected_aniline >= 2 and total_steps >= 3

    print(f"Steps with protected aniline: {steps_with_protected_aniline}/{total_steps}")
    print(f"Linear synthesis with protected aniline: {is_linear_with_protected_aniline}")

    return is_linear_with_protected_aniline
