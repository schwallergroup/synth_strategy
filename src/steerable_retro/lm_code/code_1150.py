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
    This function detects if heterocycle formation occurs in the late stage of synthesis.
    It specifically looks for benzoxazole formation in the second half of the synthesis.
    """
    benzoxazole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nco2")
    max_depth = 0
    heterocycle_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, heterocycle_formation_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if benzoxazole is in product but not in reactants
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(benzoxazole_pattern):
                benzoxazole_in_reactants = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(benzoxazole_pattern):
                        benzoxazole_in_reactants = True
                        break

                if not benzoxazole_in_reactants:
                    print(f"Benzoxazole formation detected at depth {depth}")
                    heterocycle_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if heterocycle formation occurs in the second half of synthesis
    if heterocycle_formation_depth is not None and max_depth > 0:
        relative_position = heterocycle_formation_depth / max_depth
        print(f"Heterocycle formation relative position: {relative_position:.2f}")
        return (
            relative_position <= 0.5
        )  # Second half of synthesis (remember depth increases as we go backward)
    return False
