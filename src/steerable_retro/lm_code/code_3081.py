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
    This function detects a synthetic strategy involving sequential formation of multiple sulfonamide groups.
    """
    sulfonamide_formation_count = 0
    sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#6]")

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check if sulfonamide group is formed in this reaction
                    reactants_matches = len(reactants_mol.GetSubstructMatches(sulfonamide_pattern))
                    product_matches = len(product_mol.GetSubstructMatches(sulfonamide_pattern))

                    if product_matches > reactants_matches:
                        sulfonamide_formation_count += 1
                        print(f"Detected sulfonamide formation at reaction with SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if multiple sulfonamide formations are detected
    result = sulfonamide_formation_count >= 2
    print(
        f"Sequential sulfonamide formation strategy detected: {result} (count: {sulfonamide_formation_count})"
    )
    return result
