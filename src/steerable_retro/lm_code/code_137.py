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
    This function detects a linear synthesis strategy building a complex heterocyclic system
    through sequential fragment additions.
    """
    # Track the number of sequential reactions
    reaction_count = 0
    heterocycle_reactions = 0

    def dfs_traverse(node):
        nonlocal reaction_count, heterocycle_reactions

        if node["type"] == "reaction":
            reaction_count += 1

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if product contains heterocyclic systems
                heterocycle_patterns = [
                    Chem.MolFromSmarts("[#6]1[#6][#7][#6][#6][#6]1"),  # pyridine
                    Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6]2[#6]1[#7][#6][#7]2"),  # benzimidazole
                    Chem.MolFromSmarts(
                        "[#7]1[#6][#16][#6]2[#7][#6][#7][#6][#6]21"
                    ),  # complex heterocycle
                ]

                product_mol = Chem.MolFromSmiles(product)
                if product_mol is not None:
                    for pattern in heterocycle_patterns:
                        if product_mol.HasSubstructMatch(pattern):
                            heterocycle_reactions += 1
                            break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Linear heterocycle synthesis if:
    # 1. Multiple reactions (at least 3)
    # 2. Most reactions involve heterocycles
    is_linear_heterocycle_synthesis = (
        reaction_count >= 3 and heterocycle_reactions >= reaction_count * 0.7
    )

    if is_linear_heterocycle_synthesis:
        print(
            f"Detected linear heterocycle synthesis with {reaction_count} reactions, {heterocycle_reactions} involving heterocycles"
        )

    return is_linear_heterocycle_synthesis
