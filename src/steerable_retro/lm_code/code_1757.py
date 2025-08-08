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
    Detects if the synthesis involves a sequence of functional group interconversions,
    specifically methyl → carboxylic acid → ester transformations.
    """
    # Initialize tracking variables
    methyl_to_acid = False
    acid_to_ester = False

    def dfs_traverse(node, depth=0):
        nonlocal methyl_to_acid, acid_to_ester

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for methyl to carboxylic acid transformation
            methyl_pattern = Chem.MolFromSmarts("[#6]-[#6;H3]")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8;H1]")

            has_methyl_reactant = False
            has_acid_product = False

            # Check reactants for methyl
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(methyl_pattern):
                        has_methyl_reactant = True
                except:
                    continue

            # Check product for carboxylic acid
            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    has_acid_product = True
            except:
                pass

            if has_methyl_reactant and has_acid_product:
                methyl_to_acid = True
                print(f"Detected methyl to carboxylic acid transformation at depth {depth}")

            # Check for carboxylic acid to ester transformation
            ester_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8]-[#6]")

            has_acid_reactant = False
            has_ester_product = False

            # Check reactants for carboxylic acid
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_acid_reactant = True
                except:
                    continue

            # Check product for ester
            try:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(ester_pattern):
                    has_ester_product = True
            except:
                pass

            if has_acid_reactant and has_ester_product:
                acid_to_ester = True
                print(f"Detected carboxylic acid to ester transformation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if both transformations are present
    strategy_present = methyl_to_acid and acid_to_ester
    print(f"Functional group interconversion strategy detected: {strategy_present}")
    return strategy_present
