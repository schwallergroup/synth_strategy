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
    This function detects if a nitrogen heterocycle (specifically indazole) is formed
    in the late stage of the synthesis (first half of the synthesis tree).
    """
    # Track if we found heterocycle formation
    heterocycle_formation_found = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and depth <= 3:  # Focus on late-stage reactions (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain indazole pattern
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                if reactants_mol:
                    reactants_has_indazole = reactants_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]1[#6][#6][#7][#7][#6]1")
                    )
                else:
                    reactants_has_indazole = False

                # Check if product contains indazole pattern
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    product_has_indazole = product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]1[#6][#6][#7][#7][#6]1")
                    )
                else:
                    product_has_indazole = False

                # If indazole is formed in this reaction
                if not reactants_has_indazole and product_has_indazole:
                    print(f"Detected indazole formation at depth {depth}")
                    heterocycle_formation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it late-stage if it happens in first half of synthesis
    if heterocycle_formation_found:
        print(f"Heterocycle formation found in late stage (max depth: {max_depth})")

    return heterocycle_formation_found
