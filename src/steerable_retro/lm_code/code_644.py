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
    This function detects a strategy involving heterocycle formation,
    specifically isoxazole formation from nitrile and hydroxylamine.
    """
    found_isoxazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_isoxazole_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitrile in reactants
            nitrile_present = False
            hydroxylamine_present = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    nitrile_pattern = Chem.MolFromSmarts("[#6]#[N]")
                    if mol.HasSubstructMatch(nitrile_pattern):
                        nitrile_present = True
                        print("Found nitrile in reactants")

                    hydroxylamine_pattern = Chem.MolFromSmarts("[NH2][OH]")
                    if mol.HasSubstructMatch(hydroxylamine_pattern):
                        hydroxylamine_present = True
                        print("Found hydroxylamine in reactants")

            # Check for isoxazole in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                isoxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#8][#6][#6]1")
                if product_mol.HasSubstructMatch(isoxazole_pattern):
                    # If we have both nitrile and hydroxylamine in reactants and isoxazole in product
                    if nitrile_present and hydroxylamine_present:
                        found_isoxazole_formation = True
                        print("Found isoxazole formation from nitrile and hydroxylamine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_isoxazole_formation
