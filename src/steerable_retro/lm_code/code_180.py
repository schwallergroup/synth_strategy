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
    This function detects a synthetic sequence involving ester reduction to alcohol
    followed by conversion to bromide.
    """
    # Track transformations
    ester_to_alcohol = False
    alcohol_to_bromide = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_to_alcohol, alcohol_to_bromide

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester reduction to alcohol
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
                alcohol_pattern = Chem.MolFromSmarts("[C][OH]")

                # Check for alcohol to bromide conversion
                alcohol_pattern2 = Chem.MolFromSmarts("[c][CH2][OH]")
                bromide_pattern = Chem.MolFromSmarts("[c][CH2][Br]")

                # Check reactants and products
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Ester to alcohol check
                        if (
                            reactant_mol.HasSubstructMatch(ester_pattern)
                            and product_mol.HasSubstructMatch(alcohol_pattern)
                            and not reactant_mol.HasSubstructMatch(alcohol_pattern)
                        ):
                            ester_to_alcohol = True
                            print("Found ester reduction to alcohol")

                        # Alcohol to bromide check
                        if (
                            reactant_mol.HasSubstructMatch(alcohol_pattern2)
                            and product_mol.HasSubstructMatch(bromide_pattern)
                            and not reactant_mol.HasSubstructMatch(bromide_pattern)
                        ):
                            alcohol_to_bromide = True
                            print("Found alcohol to bromide conversion")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found both transformations
    return ester_to_alcohol and alcohol_to_bromide
