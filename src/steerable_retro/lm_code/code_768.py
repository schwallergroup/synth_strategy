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
    Detects if the route follows a linear synthesis strategy without convergent branches.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are more than 2 reactants, it might be convergent
            if len(reactants_smiles) > 2:
                # Check if the additional reactants are small molecules (likely reagents)
                reactant_sizes = [
                    len(Chem.MolFromSmiles(r).GetAtoms())
                    for r in reactants_smiles
                    if Chem.MolFromSmiles(r)
                ]
                # Sort sizes in descending order
                reactant_sizes.sort(reverse=True)

                # If the third largest reactant has more than 5 atoms, consider it convergent
                if len(reactant_sizes) > 2 and reactant_sizes[2] > 5:
                    print(
                        f"Convergent synthesis detected at depth {depth} with {len(reactants_smiles)} significant reactants"
                    )
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return is_linear
