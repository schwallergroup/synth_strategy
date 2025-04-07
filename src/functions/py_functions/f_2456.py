#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects the conversion of an amide to an ester in the synthesis route.
    """
    amide_to_ester = False

    def dfs_traverse(node):
        nonlocal amide_to_ester

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide in reactants
                amide_pattern = Chem.MolFromSmarts("[#6]-[#7]-C(=O)-[#6]")
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8]")
                ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-C(=O)-[#6]")

                has_amide = False
                has_alcohol = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(amide_pattern):
                                has_amide = True
                            if mol.HasSubstructMatch(alcohol_pattern):
                                has_alcohol = True
                    except:
                        continue

                # Check if product has ester
                if has_amide and has_alcohol:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol and prod_mol.HasSubstructMatch(ester_pattern):
                            amide_to_ester = True
                            print("Found amide to ester conversion")
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_to_ester
