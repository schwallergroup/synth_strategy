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
    This function detects if the synthetic route employs a Weinreb amide formation
    followed by conversion to ketone.
    """
    # Track Weinreb amide formation and conversion
    weinreb_amide_formed = False
    weinreb_to_ketone = False

    def dfs_traverse(node):
        nonlocal weinreb_amide_formed, weinreb_to_ketone

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Weinreb amide formation
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            weinreb_amine_pattern = Chem.MolFromSmarts("[N]([CH3])[O][CH3]")
            weinreb_amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")

            # Check for Weinreb amide formation
            acid_in_reactants = False
            weinreb_in_reactants = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(carboxylic_acid_pattern):
                            acid_in_reactants = True
                        if mol.HasSubstructMatch(weinreb_amine_pattern):
                            weinreb_in_reactants = True
                except:
                    continue

            if acid_in_reactants and weinreb_in_reactants:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(weinreb_amide_pattern):
                        print("Found Weinreb amide formation")
                        weinreb_amide_formed = True
                except:
                    pass

            # Check for Weinreb amide to ketone conversion
            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[#6]")

            weinreb_in_reactants = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(weinreb_amide_pattern):
                        weinreb_in_reactants = True
                        break
                except:
                    continue

            if weinreb_in_reactants:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                        print("Found Weinreb amide to ketone conversion")
                        weinreb_to_ketone = True
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both Weinreb amide formation and conversion to ketone are found
    return weinreb_amide_formed and weinreb_to_ketone
