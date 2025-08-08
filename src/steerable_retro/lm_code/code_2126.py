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
    This function detects a synthetic strategy involving multiple C-O bond formations
    (esterifications and etherifications) in a linear synthesis.
    """
    co_bond_formation_count = 0

    def dfs_traverse(node):
        nonlocal co_bond_formation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for C-O bond formation (esterification or etherification)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for alcohol or carboxylic acid in reactants
            alcohol_pattern = Chem.MolFromSmarts("[OH]")
            acid_pattern = Chem.MolFromSmarts("[OH]C=O")

            # Look for ester or ether in product
            ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-C(=O)")
            ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

            has_alcohol_or_acid = any(
                mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(acid_pattern)
                for mol in reactant_mols
                if mol
            )

            has_ester_or_ether = product_mol and (
                product_mol.HasSubstructMatch(ester_pattern)
                or product_mol.HasSubstructMatch(ether_pattern)
            )

            if has_alcohol_or_acid and has_ester_or_ether:
                co_bond_formation_count += 1
                print(f"Detected C-O bond formation in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if there are at least 2 C-O bond formations
    return co_bond_formation_count >= 2
