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
    This function detects a synthesis strategy that involves sequential functional group
    transformations at a benzylic carbon position (CH3 → CH2Br → CH2CN → COOH).
    """
    # Track if we've found the required functional group transformations
    found_benzylic_methyl = False
    found_benzyl_bromide = False
    found_benzyl_nitrile = False
    found_benzylic_carboxylic_acid = False

    # SMARTS patterns for the functional groups
    benzylic_methyl_pattern = Chem.MolFromSmarts("[c]-[CH3]")
    benzyl_bromide_pattern = Chem.MolFromSmarts("[c]-[CH2][Br]")
    benzyl_nitrile_pattern = Chem.MolFromSmarts("[c]-[CH2][C]#[N]")
    benzylic_carboxylic_acid_pattern = Chem.MolFromSmarts("[c]-[CH2][C](=[O])[OH]")

    def dfs_traverse(node):
        nonlocal found_benzylic_methyl, found_benzyl_bromide, found_benzyl_nitrile, found_benzylic_carboxylic_acid

        if node["type"] == "mol":
            # Check for functional groups in molecule nodes
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(benzylic_methyl_pattern):
                    found_benzylic_methyl = True
                    print(f"Found benzylic methyl in: {node['smiles']}")
                if mol.HasSubstructMatch(benzyl_bromide_pattern):
                    found_benzyl_bromide = True
                    print(f"Found benzyl bromide in: {node['smiles']}")
                if mol.HasSubstructMatch(benzyl_nitrile_pattern):
                    found_benzyl_nitrile = True
                    print(f"Found benzyl nitrile in: {node['smiles']}")
                if mol.HasSubstructMatch(benzylic_carboxylic_acid_pattern):
                    found_benzylic_carboxylic_acid = True
                    print(f"Found benzylic carboxylic acid in: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found at least 3 of the 4 functional groups
    # This indicates the benzylic carbon functionalization strategy
    functional_groups_found = sum(
        [
            found_benzylic_methyl,
            found_benzyl_bromide,
            found_benzyl_nitrile,
            found_benzylic_carboxylic_acid,
        ]
    )

    strategy_detected = functional_groups_found >= 3
    print(f"Benzylic carbon functionalization strategy detected: {strategy_detected}")
    return strategy_detected
