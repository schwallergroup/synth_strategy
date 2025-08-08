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
    This function detects a sequence of functional group interconversions
    (e.g., methyl → hydroxymethyl → aldehyde).
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # For simplicity, consider only the first reactant
                    reactant_mol = Chem.MolFromSmiles(reactants[0])
                    product_mol = Chem.MolFromSmiles(product)

                    # Define patterns for functional groups
                    methyl_pattern = Chem.MolFromSmarts("[#6][CH3]")
                    hydroxymethyl_pattern = Chem.MolFromSmarts("[#6][CH2][OH]")
                    aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")

                    # Check for transformations
                    if reactant_mol.HasSubstructMatch(
                        methyl_pattern
                    ) and product_mol.HasSubstructMatch(hydroxymethyl_pattern):
                        transformations.append("methyl_to_hydroxymethyl")
                        print("Detected: methyl → hydroxymethyl")

                    if reactant_mol.HasSubstructMatch(
                        hydroxymethyl_pattern
                    ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                        transformations.append("hydroxymethyl_to_aldehyde")
                        print("Detected: hydroxymethyl → aldehyde")

                    if reactant_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ) and product_mol.HasSubstructMatch(hydroxymethyl_pattern):
                        transformations.append("aldehyde_to_hydroxymethyl")
                        print("Detected: aldehyde → hydroxymethyl")
                except:
                    print("Error processing molecules for FG transformation detection")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a sequence of transformations
    return len(transformations) >= 2
