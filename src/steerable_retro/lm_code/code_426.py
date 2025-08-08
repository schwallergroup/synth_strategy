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
    Detects if the synthesis involves a multicomponent reaction,
    particularly focusing on reactions with 3+ reactants forming complex products.
    """
    has_multicomponent_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_multicomponent_reaction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if there are 3 or more distinct reactants
                if len(reactants) >= 3:
                    # Verify these are truly distinct molecules, not just fragments
                    distinct_reactants = set()
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # Use a simplified representation to avoid counting very similar molecules
                                # as distinct (e.g., different salts of the same compound)
                                simplified = Chem.MolToSmiles(mol, isomericSmiles=False)
                                distinct_reactants.add(simplified)
                        except:
                            continue

                    if len(distinct_reactants) >= 3:
                        has_multicomponent_reaction = True
                        print(
                            f"Detected multicomponent reaction with {len(distinct_reactants)} distinct reactants at depth {depth}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_multicomponent_reaction
