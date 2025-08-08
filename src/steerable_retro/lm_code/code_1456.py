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
    This function detects a strategy involving sequential transformation:
    isothiocyanate formation followed by thiourea formation.
    """
    isothiocyanate_formed = False
    thiourea_from_isothiocyanate = False

    def dfs_traverse(node):
        nonlocal isothiocyanate_formed, thiourea_from_isothiocyanate

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for isothiocyanate formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    isothiocyanate_pattern = Chem.MolFromSmarts("[#6]-[N]=[C]=[S]")
                    if product_mol.HasSubstructMatch(isothiocyanate_pattern):
                        isothiocyanate_formed = True
                        print("Found isothiocyanate formation")

                # Check for thiourea formation from isothiocyanate
                thiourea_pattern = Chem.MolFromSmarts("[#6]-[NH]-[C](=[S])-[NH2]")
                if product_mol and product_mol.HasSubstructMatch(thiourea_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(isothiocyanate_pattern):
                                thiourea_from_isothiocyanate = True
                                print("Found thiourea formation from isothiocyanate")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return isothiocyanate_formed and thiourea_from_isothiocyanate
