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
    This function detects if the synthetic route employs Fischer indole synthesis
    to construct an indole ring system.
    """
    fischer_indole_detected = False

    def dfs_traverse(node):
        nonlocal fischer_indole_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain hydrazine derivative
            hydrazine_pattern = Chem.MolFromSmarts("[N]-[N]")
            carbonyl_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]")
            indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")

            # Check reactants for hydrazine and carbonyl patterns
            hydrazine_found = False
            carbonyl_found = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(hydrazine_pattern):
                        hydrazine_found = True
                    if mol and mol.HasSubstructMatch(carbonyl_pattern):
                        carbonyl_found = True
                except:
                    continue

            # Check product for indole pattern
            try:
                product_mol = Chem.MolFromSmiles(product)
                indole_found = product_mol and product_mol.HasSubstructMatch(indole_pattern)
            except:
                indole_found = False

            # If reactants contain hydrazine and carbonyl, and product contains indole, it's likely Fischer indole synthesis
            if hydrazine_found and carbonyl_found and indole_found:
                print("Fischer indole synthesis detected")
                fischer_indole_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fischer_indole_detected
