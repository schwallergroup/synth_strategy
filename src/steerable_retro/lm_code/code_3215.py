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
    This function detects a strategy where halogen substituents (Cl, F) are preserved
    throughout a linear functionalization sequence.
    """
    # Initialize tracking variables
    depths_with_reactions = []
    all_reactions_preserve_halogens = True

    # Define SMARTS patterns for halogens
    chloro_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    fluoro_pattern = Chem.MolFromSmarts("[#6]-[F]")

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_reactions, all_reactions_preserve_halogens

        if node["type"] == "reaction":
            depths_with_reactions.append(depth)

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has both Cl and F
            product_mol = Chem.MolFromSmiles(product_smiles)
            if not (
                product_mol
                and product_mol.HasSubstructMatch(chloro_pattern)
                and product_mol.HasSubstructMatch(fluoro_pattern)
            ):
                all_reactions_preserve_halogens = False
                print(f"Reaction at depth {depth} does not preserve both halogens")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if there are multiple reactions and all preserve halogens
    strategy_present = len(depths_with_reactions) >= 2 and all_reactions_preserve_halogens

    print(f"Halogen preserving functionalization strategy detected: {strategy_present}")
    return strategy_present
