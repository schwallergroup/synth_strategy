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
    This function detects a synthesis strategy involving the transformation
    of a nitrile to an amide via a carboxylic acid intermediate.
    """
    # Track functional group transformations
    has_nitrile_to_acid = False
    has_acid_to_amide = False

    def dfs_traverse(node):
        nonlocal has_nitrile_to_acid, has_acid_to_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define patterns
                nitrile_patt = Chem.MolFromSmarts("[#6]-[#6]#[#7]")
                carboxylic_patt = Chem.MolFromSmarts("[#6]-[#6](=O)-[#8;H1]")
                amide_patt = Chem.MolFromSmarts("[#6]-[#6](=O)-[#7]")

                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Check for nitrile to carboxylic acid transformation
                    if reactant_mol.HasSubstructMatch(
                        nitrile_patt
                    ) and product_mol.HasSubstructMatch(carboxylic_patt):
                        has_nitrile_to_acid = True
                        print("Found nitrile to carboxylic acid transformation")

                    # Check for carboxylic acid to amide transformation
                    if reactant_mol.HasSubstructMatch(
                        carboxylic_patt
                    ) and product_mol.HasSubstructMatch(amide_patt):
                        has_acid_to_amide = True
                        print("Found carboxylic acid to amide transformation")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = has_nitrile_to_acid and has_acid_to_amide
    print(f"Nitrile to amide via acid strategy detected: {result}")
    return result
