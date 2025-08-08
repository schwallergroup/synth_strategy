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
    This function detects a synthesis strategy involving bidirectional
    transformations between alcohol and aldehyde functional groups.
    """
    # Track oxidation and reduction steps
    has_alcohol_to_aldehyde = False
    has_aldehyde_to_alcohol = False

    def dfs_traverse(node):
        nonlocal has_alcohol_to_aldehyde, has_aldehyde_to_alcohol

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol to aldehyde transformation (oxidation)
                alcohol_patt = Chem.MolFromSmarts("[#6]-[#8;H1]")
                aldehyde_patt = Chem.MolFromSmarts("[#6;H1]=O")

                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    if reactant_mol.HasSubstructMatch(
                        alcohol_patt
                    ) and product_mol.HasSubstructMatch(aldehyde_patt):
                        has_alcohol_to_aldehyde = True
                        print("Found alcohol to aldehyde transformation")

                    if reactant_mol.HasSubstructMatch(
                        aldehyde_patt
                    ) and product_mol.HasSubstructMatch(alcohol_patt):
                        has_aldehyde_to_alcohol = True
                        print("Found aldehyde to alcohol transformation")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = has_alcohol_to_aldehyde and has_aldehyde_to_alcohol
    print(f"Bidirectional alcohol-aldehyde strategy detected: {result}")
    return result
