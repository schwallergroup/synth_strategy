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
    This function detects if the synthesis involves a sequence of functional group
    interconversions: ester -> acid -> amide or similar patterns.
    """
    # Track functional group transformations
    transformations = []
    result = False

    def dfs_traverse(node):
        nonlocal transformations, result

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Define patterns for functional groups
                ester_pattern = Chem.MolFromSmarts("[C$(C=O)][O][C]")
                acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")
                amide_pattern = Chem.MolFromSmarts("[C$(C=O)][N]")

                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                product_mol = Chem.MolFromSmiles(product_smiles)

                if reactants_mol and product_mol:
                    # Check for ester hydrolysis
                    if reactants_mol.HasSubstructMatch(
                        ester_pattern
                    ) and product_mol.HasSubstructMatch(acid_pattern):
                        transformations.append("ester_to_acid")
                        print("Ester to acid transformation detected")

                    # Check for esterification
                    if reactants_mol.HasSubstructMatch(
                        acid_pattern
                    ) and product_mol.HasSubstructMatch(ester_pattern):
                        transformations.append("acid_to_ester")
                        print("Acid to ester transformation detected")

                    # Check for amide formation
                    if reactants_mol.HasSubstructMatch(
                        acid_pattern
                    ) and product_mol.HasSubstructMatch(amide_pattern):
                        transformations.append("acid_to_amide")
                        print("Acid to amide transformation detected")
            except Exception as e:
                print(f"Error in functional group detection: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check for specific sequences
    if "ester_to_acid" in transformations and "acid_to_amide" in transformations:
        print("Detected ester -> acid -> amide sequence")
        result = True
    if "acid_to_ester" in transformations and "ester_to_acid" in transformations:
        print("Detected acid -> ester -> acid sequence")
        result = True

    return result
