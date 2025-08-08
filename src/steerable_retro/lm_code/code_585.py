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
    This function detects conversion of alkyl halide to alcohol in the early stage of synthesis.
    It looks for reactions where a bromomethyl ketone is converted to a hydroxymethyl ketone.
    """
    halide_to_alcohol_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal halide_to_alcohol_detected

        if node["type"] == "reaction" and depth >= 2:  # Early stage (far from final product)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Define patterns
            bromomethyl_ketone = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]-[#35]")
            hydroxymethyl_ketone = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]-[#8]")

            # Check if bromomethyl ketone is in reactants
            bromo_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(bromomethyl_ketone):
                    bromo_in_reactants = True
                    break

            # Check if hydroxymethyl ketone is in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            hydroxy_in_product = product_mol and product_mol.HasSubstructMatch(hydroxymethyl_ketone)

            if bromo_in_reactants and hydroxy_in_product:
                print(f"Halide to alcohol conversion detected at depth {depth}")
                halide_to_alcohol_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return halide_to_alcohol_detected
