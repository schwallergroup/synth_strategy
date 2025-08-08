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
    This function detects a synthesis strategy that utilizes phenylacetylene as a key building block,
    particularly through Sonogashira coupling.
    """
    # Initialize tracking variables
    has_phenylacetylene = False
    has_sonogashira = False

    def dfs_traverse(node):
        nonlocal has_phenylacetylene, has_sonogashira

        if node["type"] == "mol":
            # Check if this molecule is phenylacetylene
            smiles = node.get("smiles", "")
            if smiles and node.get("in_stock", False):
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    phenylacetylene_pattern = Chem.MolFromSmarts("c1ccccc1C#C")
                    if mol.HasSubstructMatch(phenylacetylene_pattern):
                        has_phenylacetylene = True
                        print("Detected phenylacetylene as starting material")

        elif node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Sonogashira coupling
            if "#" in product_smiles:
                # Check if one reactant has halogen and another has terminal alkyne
                has_halogen = any(re.search(r"Br|I|Cl", r) for r in reactants_smiles)
                has_alkyne = any("#" in r for r in reactants_smiles)

                if has_halogen and has_alkyne:
                    has_sonogashira = True
                    print("Detected Sonogashira coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Check if both required elements are present
    strategy_detected = has_phenylacetylene and has_sonogashira

    return strategy_detected
