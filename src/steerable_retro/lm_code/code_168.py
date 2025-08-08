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
    Detects a strategy involving nitro-substituted aromatic compounds throughout the synthesis.
    """
    depths_with_nitro = set()

    # SMARTS pattern for nitro group
    nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(nitro_pattern):
                depths_with_nitro.add(depth)
                print(f"Found nitro group at depth {depth}")

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            for smiles in reactants_smiles + [product_smiles]:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(nitro_pattern):
                    depths_with_nitro.add(depth)
                    print(f"Found nitro group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if nitro groups are found at multiple depths
    strategy_present = len(depths_with_nitro) >= 2

    if strategy_present:
        print(f"Detected nitro-containing aromatic strategy across depths: {depths_with_nitro}")

    return strategy_present
