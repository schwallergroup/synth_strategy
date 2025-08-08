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
    This function detects the conversion of alcohols to chlorides in the synthetic route.
    """
    conversion_detected = False

    def dfs_traverse(node):
        nonlocal conversion_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain alcohol and product contains chloride
                alcohol_pattern = Chem.MolFromSmarts("[#8H1]-[#6]")
                chloride_pattern = Chem.MolFromSmarts("[Cl]-[#6]")

                reactant_has_alcohol = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(alcohol_pattern):
                        reactant_has_alcohol = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                product_has_chloride = product_mol and product_mol.HasSubstructMatch(
                    chloride_pattern
                )

                if reactant_has_alcohol and product_has_chloride:
                    conversion_detected = True
                    print("Alcohol to chloride conversion detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return conversion_detected
