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
    Detects a linear synthesis strategy that maintains a chloro-substituted aromatic ring.
    """
    chloro_aromatic_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal chloro_aromatic_count, max_reactants_per_step

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = reactants_smiles.split(".")
                num_reactants = len(reactants)
                max_reactants_per_step = max(max_reactants_per_step, num_reactants)

                product = Chem.MolFromSmiles(product_smiles)

                if product:
                    # Check for chloro-aromatic pattern
                    chloro_aromatic_pattern = Chem.MolFromSmarts("[Cl][c]")
                    if product.HasSubstructMatch(chloro_aromatic_pattern):
                        chloro_aromatic_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis typically has 1-2 reactants per step
    is_linear = max_reactants_per_step <= 2
    maintains_chloro_aromatic = chloro_aromatic_count >= 2  # Present in multiple steps

    print(f"Linear synthesis: {is_linear}, Chloro-aromatic count: {chloro_aromatic_count}")
    return is_linear and maintains_chloro_aromatic
