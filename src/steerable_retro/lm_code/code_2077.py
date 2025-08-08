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
    This function detects a strategy involving functionalization of a piperazine scaffold
    that is maintained throughout the synthesis.
    """
    # Track if we found piperazine scaffold and its functionalization
    found_piperazine = False
    found_functionalization = False

    def dfs_traverse(node):
        nonlocal found_piperazine, found_functionalization

        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule contains piperazine scaffold
            piperazine_pattern = Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#7]-[#6]-[#6]-1")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(piperazine_pattern):
                found_piperazine = True
                print(f"Found piperazine scaffold in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for functionalization reactions on piperazine
            piperazine_pattern = Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#7]-[#6]-[#6]-1")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(piperazine_pattern):
                # Check if this is a functionalization reaction
                # (We're looking for reactions that modify the piperazine but keep its core intact)
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(piperazine_pattern):
                        # If both reactant and product have piperazine, it's a functionalization
                        found_functionalization = True
                        print(f"Found piperazine functionalization in reaction: {rsmi}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found piperazine scaffold and its functionalization
    strategy_present = found_piperazine and found_functionalization

    if strategy_present:
        print("Detected strategy: Piperazine scaffold functionalization")

    return strategy_present
