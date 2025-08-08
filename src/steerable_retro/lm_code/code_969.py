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
    This function detects if the synthetic route includes borylation of an aryl bromide.
    """
    borylation_found = False

    def dfs_traverse(node):
        nonlocal borylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl bromide in reactants
            aryl_bromide_pattern = Chem.MolFromSmarts("[c][Br]")
            # Check for boronic acid in product
            boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")

            aryl_bromide_present = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_bromide_pattern):
                        aryl_bromide_present = True
                        break
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                boronic_acid_present = product_mol and product_mol.HasSubstructMatch(
                    boronic_acid_pattern
                )
            except:
                boronic_acid_present = False

            if aryl_bromide_present and boronic_acid_present:
                print("Borylation of aryl bromide detected")
                borylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return borylation_found
