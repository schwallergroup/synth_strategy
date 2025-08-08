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
    This function detects a synthetic strategy involving acid chloride transformations.
    """
    acid_chloride_transformation_detected = False

    def dfs_traverse(node):
        nonlocal acid_chloride_transformation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for acid chloride and carboxylic acid
                acid_chloride_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#17]")
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8]")

                # Check reactants and products
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                # Check for acid chloride to carboxylic acid transformation
                if product_mol and product_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    if any(
                        mol and mol.HasSubstructMatch(acid_chloride_pattern)
                        for mol in reactant_mols
                    ):
                        print("Detected acid chloride to carboxylic acid transformation")
                        acid_chloride_transformation_detected = True

                # Also check reverse (retrosynthetic direction)
                if product_mol and product_mol.HasSubstructMatch(acid_chloride_pattern):
                    if any(
                        mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                        for mol in reactant_mols
                    ):
                        print("Detected carboxylic acid to acid chloride transformation")
                        acid_chloride_transformation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return acid_chloride_transformation_detected
