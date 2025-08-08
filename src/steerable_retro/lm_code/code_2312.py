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
    This function detects a linear synthesis strategy involving
    an acyl hydrazide intermediate formation.
    """
    # Initialize tracking variables
    has_acyl_hydrazide = False
    is_linear_synthesis = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_acyl_hydrazide, is_linear_synthesis, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a linear step (one product)
            if "." in product_smiles:
                is_linear_synthesis = False
                print(f"Non-linear step detected at depth {depth}: multiple products")

            # Convert to RDKit molecules
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product:
                # Check for acyl hydrazide formation
                acyl_hydrazide_pattern = Chem.MolFromSmarts("[C](=[O])[N][N]")
                if product.HasSubstructMatch(acyl_hydrazide_pattern):
                    has_acyl_hydrazide = True
                    print(f"Detected acyl hydrazide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if conditions are met
    return is_linear_synthesis and has_acyl_hydrazide and reaction_count >= 3
