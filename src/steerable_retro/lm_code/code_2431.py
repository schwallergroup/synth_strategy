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
    Detects if the synthesis route involves multiple steps of alcohol manipulations
    (protection, functionalization, etc.)
    """
    alcohol_steps = 0

    def dfs_traverse(node):
        nonlocal alcohol_steps

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for alcohol-related transformations
            alcohol_pattern = Chem.MolFromSmarts("[O;!H0]")
            silyl_ether_pattern = Chem.MolFromSmarts("[O]-[Si]")
            ester_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
            thiocarbamate_pattern = Chem.MolFromSmarts("[O]-[C](=[S])-[O]")

            product_mol = Chem.MolFromSmiles(product_smiles)

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)

                if reactant_mol and product_mol:
                    # Check for various alcohol transformations
                    if reactant_mol.HasSubstructMatch(alcohol_pattern) and (
                        product_mol.HasSubstructMatch(silyl_ether_pattern)
                        or product_mol.HasSubstructMatch(thiocarbamate_pattern)
                    ):
                        alcohol_steps += 1
                        print("Alcohol protection or functionalization detected")

                    if reactant_mol.HasSubstructMatch(
                        ester_pattern
                    ) and product_mol.HasSubstructMatch(alcohol_pattern):
                        alcohol_steps += 1
                        print("Ester to alcohol transformation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Multiple alcohol manipulations (at least 2)
    multiple_manipulations = alcohol_steps >= 2

    if multiple_manipulations:
        print(f"Multiple alcohol manipulations detected: {alcohol_steps} steps")

    return multiple_manipulations
