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
    Detects if the synthesis includes the sequence: COOH → COOCH3 → CH2OH → CH2Cl
    """
    # Track which transformations we've seen
    transformations = {
        "acid_to_ester": False,
        "ester_to_alcohol": False,
        "alcohol_to_chloride": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid to ester transformation
            acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")

            # Check for ester to alcohol transformation
            ester_to_alcohol_pattern_reactant = Chem.MolFromSmarts("[C](=O)[O][C]")
            alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")

            # Check for alcohol to chloride transformation
            alcohol_to_chloride_pattern_reactant = Chem.MolFromSmarts("[CH2][OH]")
            chloride_pattern = Chem.MolFromSmarts("[CH2][Cl]")

            # Check acid to ester
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(acid_pattern):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(ester_pattern):
                        transformations["acid_to_ester"] = True
                        print("Found acid to ester transformation")

            # Check ester to alcohol
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(ester_to_alcohol_pattern_reactant):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(alcohol_pattern):
                        transformations["ester_to_alcohol"] = True
                        print("Found ester to alcohol transformation")

            # Check alcohol to chloride
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(alcohol_to_chloride_pattern_reactant):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(chloride_pattern):
                        transformations["alcohol_to_chloride"] = True
                        print("Found alcohol to chloride transformation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we've seen all three transformations
    return (
        transformations["acid_to_ester"]
        and transformations["ester_to_alcohol"]
        and transformations["alcohol_to_chloride"]
    )
