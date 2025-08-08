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
    Detects a synthetic strategy involving the formation of a heterocyclic ring
    in the final step of the synthesis.
    """
    # Flag to track if heterocycle formation occurs in the final step
    final_step_heterocycle_formation = False

    # Common heterocycle SMARTS patterns
    triazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#6]1")

    def dfs_traverse(node):
        nonlocal final_step_heterocycle_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is the final step (depth 0)
            if node.get("depth", 0) == 0:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check reactants and products for heterocycle patterns
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and reactants_mol:
                    # Check if product has heterocycle but reactants don't
                    if product_mol.HasSubstructMatch(
                        triazole_pattern
                    ) and not reactants_mol.HasSubstructMatch(triazole_pattern):
                        final_step_heterocycle_formation = True
                        print("Detected heterocycle formation in final step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if final_step_heterocycle_formation:
        print("Confirmed heterocycle formation strategy in final step")
    else:
        print("No heterocycle formation detected in final step")

    return final_step_heterocycle_formation
