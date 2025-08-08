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
    Detects a borylation followed by Suzuki coupling sequence.
    """
    # Initialize flags
    found_borylation = False
    found_suzuki = False

    # Track depths
    borylation_depth = -1
    suzuki_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_borylation, found_suzuki, borylation_depth, suzuki_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product or not all(reactants):
                return

            # Check for borylation
            br_pattern = Chem.MolFromSmarts("[c]-[Br]")
            boron_pattern = Chem.MolFromSmarts("[c]-[B]")

            if any(
                r.HasSubstructMatch(br_pattern) for r in reactants
            ) and product.HasSubstructMatch(boron_pattern):
                found_borylation = True
                borylation_depth = depth
                print(f"Detected borylation at depth {depth}")

            # Check for Suzuki coupling
            if any(r.HasSubstructMatch(boron_pattern) for r in reactants) and any(
                r.HasSubstructMatch(br_pattern) for r in reactants
            ):
                # Check if product has a biaryl bond
                biaryl_pattern = Chem.MolFromSmarts("c-c")
                if product.HasSubstructMatch(biaryl_pattern):
                    found_suzuki = True
                    suzuki_depth = depth
                    print(f"Detected Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both transformations were found in the correct order
    correct_sequence = found_borylation and found_suzuki and borylation_depth > suzuki_depth

    print(f"Cross-coupling sequence detected: {correct_sequence}")
    return correct_sequence
