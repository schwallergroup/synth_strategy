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
    Detects if the synthesis includes nitrogen protection (specifically Boc protection)
    in the middle of the synthesis route.
    """
    has_n_protection = False
    total_depth = 0
    protection_depth = -1

    # First pass to determine total depth
    def get_max_depth(node, current_depth=0):
        nonlocal total_depth
        total_depth = max(total_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    # Second pass to check for protection
    def dfs_traverse(node, depth=0):
        nonlocal has_n_protection, protection_depth

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Create RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not reactant_mols or not product_mol:
                return

            # Check for N-Boc protection
            nh_pattern = Chem.MolFromSmarts("[NH]")
            nboc_pattern = Chem.MolFromSmarts("[N][C](=[O])[O][C]([CH3])([CH3])[CH3]")

            if any(
                mol.HasSubstructMatch(nh_pattern) for mol in reactant_mols
            ) and product_mol.HasSubstructMatch(nboc_pattern):
                has_n_protection = True
                protection_depth = depth
                print(f"Found N-Boc protection at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Run both passes
    get_max_depth(route)
    dfs_traverse(route)

    # Check if protection is in the middle third of the synthesis
    is_mid_synthesis = False
    if has_n_protection and total_depth > 0:
        # Consider "middle" as between 25% and 75% of the total depth
        lower_bound = total_depth * 0.25
        upper_bound = total_depth * 0.75
        is_mid_synthesis = lower_bound <= protection_depth <= upper_bound
        print(f"Protection depth: {protection_depth}, Total depth: {total_depth}")
        print(f"Middle synthesis bounds: {lower_bound} to {upper_bound}")

    result = has_n_protection and is_mid_synthesis
    print(f"Mid-synthesis nitrogen protection detection result: {result}")
    return result
