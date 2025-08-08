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
    Detects a chain of at least 3 consecutive functional group interconversions
    (e.g., aldehyde → alcohol → acetate → chloride).
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0, parent_smiles=None):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Identify functional groups in the product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    fg_type = identify_functional_group(product_mol)
                    if fg_type:
                        transformations.append((depth, fg_type))
                        print(f"Found functional group {fg_type} at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, node.get("smiles"))

    def identify_functional_group(mol):
        """Identify key functional groups in the molecule"""
        # Aldehyde
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][#6](=[O])[H]")):
            return "aldehyde"
        # Alcohol
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][#6][OH]")):
            return "alcohol"
        # Acetate
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][#6][O][C](=[O])[#6]")):
            return "acetate"
        # Chloride
        elif mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][#6][Cl]")):
            return "chloride"
        return None

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth
    transformations.sort(key=lambda x: x[0])

    # Check for consecutive transformations
    consecutive_count = 1
    max_consecutive = 1
    for i in range(1, len(transformations)):
        if transformations[i][1] != transformations[i - 1][1]:
            consecutive_count += 1
            max_consecutive = max(max_consecutive, consecutive_count)
        else:
            consecutive_count = 1

    print(f"Maximum consecutive functional group transformations: {max_consecutive}")
    return max_consecutive >= 3
