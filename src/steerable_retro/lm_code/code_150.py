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
    Detects if the synthesis follows a linear strategy (no convergent steps).

    A linear synthesis has sequential steps where each reaction has only one
    significant reactant (the previous intermediate). Convergent synthesis
    involves steps where multiple complex reactants are combined.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Filter out small molecules and common reagents
                # Consider a reactant significant if it has more than 6 atoms
                significant_reactants = []
                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol and len(mol.GetAtoms()) > 6:
                        significant_reactants.append(r)

                if len(significant_reactants) > 1:
                    print(
                        f"Found convergent step with {len(significant_reactants)} significant reactants"
                    )
                    for i, r in enumerate(significant_reactants):
                        print(f"  Reactant {i+1}: {r}")
                    is_linear = False
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return is_linear
