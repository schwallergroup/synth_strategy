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
    This function detects a synthetic strategy involving multiple steps where
    alcohols are functionalized, particularly through etherification reactions.
    """
    alcohol_functionalization_count = 0

    def dfs_traverse(node):
        nonlocal alcohol_functionalization_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for alcohol functionalization
            alcohol_pattern = Chem.MolFromSmarts("[OH]")
            ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

            has_alcohol_reactant = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(alcohol_pattern):
                        has_alcohol_reactant = True
                        break
                except:
                    continue

            if has_alcohol_reactant:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # Check if alcohol is converted to ether
                        if product_mol.HasSubstructMatch(ether_pattern):
                            alcohol_functionalization_count += 1
                            print(
                                f"Detected alcohol functionalization in reaction {node.get('metadata', {}).get('ID', '')}"
                            )
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if multiple alcohol functionalizations are detected
    return alcohol_functionalization_count >= 2
