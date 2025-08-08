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
    This function detects a linear synthesis strategy where fragments are added
    sequentially rather than in a convergent manner.
    """
    # Track the number of reactions with multiple fragments
    multi_fragment_reactions = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal multi_fragment_reactions, total_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]

                # Count the number of reactants
                reactant_count = len(reactants_smiles.split("."))

                if reactant_count >= 2:
                    multi_fragment_reactions += 1
                    print(f"Detected multi-fragment reaction with {reactant_count} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If most reactions involve multiple fragments and we have a reasonable number of reactions,
    # consider it a linear fragment assembly strategy
    return total_reactions >= 3 and multi_fragment_reactions >= 3
