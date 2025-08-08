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
    This function detects if the synthesis involves multiple nitrogen oxidation state changes
    (nitro → amine → amide sequence).
    """
    # Define patterns for different nitrogen oxidation states
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    amine_pattern = Chem.MolFromSmarts("[NH2]")
    amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")

    # Track if we've seen each state
    nitro_seen = False
    amine_seen = False
    amide_seen = False

    def dfs_traverse(node):
        nonlocal nitro_seen, amine_seen, amide_seen

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(nitro_pattern):
                    nitro_seen = True
                    print(f"Nitro group found in: {node['smiles']}")
                if mol.HasSubstructMatch(amine_pattern):
                    amine_seen = True
                    print(f"Amine group found in: {node['smiles']}")
                if mol.HasSubstructMatch(amide_pattern):
                    amide_seen = True
                    print(f"Amide group found in: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have the sequence nitro → amine → amide
    has_oxidation_state_changes = nitro_seen and amine_seen and amide_seen
    print(f"Nitrogen oxidation state changes strategy detected: {has_oxidation_state_changes}")
    return has_oxidation_state_changes
