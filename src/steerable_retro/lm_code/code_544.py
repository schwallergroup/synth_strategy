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
    Detects a functional group interconversion sequence: ketone → alcohol → ketone.
    """
    ketone_to_alcohol = False
    alcohol_to_ketone = False

    def dfs_traverse(node, depth=0):
        nonlocal ketone_to_alcohol, alcohol_to_ketone

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants[0])
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                # Pattern for ketone
                ketone_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[#6]")
                # Pattern for secondary alcohol
                alcohol_pattern = Chem.MolFromSmarts("[#6][#6]([#8H])[#6]")

                # Check for ketone → alcohol (reduction)
                if (
                    reactant_mol.HasSubstructMatch(ketone_pattern)
                    and product_mol.HasSubstructMatch(alcohol_pattern)
                    and not product_mol.HasSubstructMatch(ketone_pattern)
                ):
                    ketone_to_alcohol = True
                    print("Found ketone to alcohol reduction")

                # Check for alcohol → ketone (oxidation)
                if (
                    reactant_mol.HasSubstructMatch(alcohol_pattern)
                    and product_mol.HasSubstructMatch(ketone_pattern)
                    and not reactant_mol.HasSubstructMatch(ketone_pattern)
                ):
                    alcohol_to_ketone = True
                    print("Found alcohol to ketone oxidation")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ketone_to_alcohol and alcohol_to_ketone
