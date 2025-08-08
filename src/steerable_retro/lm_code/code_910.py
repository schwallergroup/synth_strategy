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
    Detects a synthetic strategy that maintains a stereocenter throughout multiple transformations.
    """
    # Track reactions with stereocenter
    reactions_with_stereocenter = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal reactions_with_stereocenter, total_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                try:
                    # Check for presence of stereocenter in product
                    p_mol = Chem.MolFromSmiles(product)
                    if p_mol:
                        # Look for @ or @@ in the SMILES string to indicate stereochemistry
                        if "@" in product:
                            reactions_with_stereocenter += 1
                            print(f"Detected stereocenter in reaction product: {product}")

                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if stereocenter is preserved in at least 70% of reactions and there are at least 3 reactions
    return total_reactions >= 3 and (reactions_with_stereocenter / total_reactions >= 0.7)
