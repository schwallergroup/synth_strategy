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
    This function detects a sequence of oxidation steps leading to amide formation:
    aldehyde -> carboxylic acid -> amide
    """
    # Track if we've seen each step in the sequence
    seen_aldehyde = False
    seen_carboxylic_acid = False
    seen_amide = False

    def dfs_traverse(node):
        nonlocal seen_aldehyde, seen_carboxylic_acid, seen_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Check for functional groups
                        aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                        amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")

                        if product_mol.HasSubstructMatch(aldehyde_pattern):
                            seen_aldehyde = True
                        if product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                            seen_carboxylic_acid = True
                        if product_mol.HasSubstructMatch(amide_pattern):
                            seen_amide = True
                except:
                    print("Error processing molecule in oxidation_sequence_to_amide")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if we've seen the complete sequence
    if seen_aldehyde and seen_carboxylic_acid and seen_amide:
        print("Detected oxidation sequence: aldehyde -> carboxylic acid -> amide")
        return True
    return False
