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
    Detects a convergent synthesis strategy where multiple complex fragments
    are joined together rather than linear elaboration of a single core.
    """
    # Track fragment couplings at different depths
    fragment_couplings = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            reactants = reactants_part.split(".")

            # Check if this is a fragment coupling (multiple complex reactants)
            if len(reactants) >= 2:
                complex_reactants = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).GetNumAtoms() > 8
                )

                if complex_reactants >= 2:
                    fragment_couplings.append((depth, complex_reactants))
                    print(
                        f"Found fragment coupling at depth {depth} with {complex_reactants} complex reactants"
                    )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze the pattern of fragment couplings
    # Convergent synthesis typically has multiple fragment couplings at different depths
    is_convergent = len(fragment_couplings) >= 2 and len(set(d for d, _ in fragment_couplings)) >= 2

    print(f"Convergent synthesis with multiple fragment couplings detected: {is_convergent}")
    return is_convergent
