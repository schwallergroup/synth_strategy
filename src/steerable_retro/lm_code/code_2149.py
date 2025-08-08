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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If more than 2 reactants, it's potentially a convergent step
            if len(reactants_smiles) > 2:
                # Check if they are actual distinct fragments (not just reagents)
                significant_fragments = 0
                for reactant in reactants_smiles:
                    mol = Chem.MolFromSmiles(reactant)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 3
                    ):  # Arbitrary threshold to exclude small reagents
                        significant_fragments += 1

                if significant_fragments > 2:
                    print(
                        f"Convergent step detected at depth {depth} with {significant_fragments} significant fragments"
                    )
                    is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear
