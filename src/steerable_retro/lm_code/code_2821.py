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
    This function detects phenol protection (as benzyl ether) followed by deprotection.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for protection (phenol + benzyl source -> benzyl protected phenol)
            reactants = reactants_part.split(".")
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol:
                benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")
                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                    # Check if reactants include phenol
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                            if mol.HasSubstructMatch(phenol_pattern):
                                protection_found = True
                                print("Found phenol protection with benzyl group")

            # Check for deprotection (benzyl protected phenol -> phenol)
            reactant_mol = Chem.MolFromSmiles(reactants_part)
            if reactant_mol:
                benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")
                if reactant_mol.HasSubstructMatch(benzyl_ether_pattern):
                    # Check if product includes phenol
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol:
                        phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                        if product_mol.HasSubstructMatch(phenol_pattern):
                            deprotection_found = True
                            print("Found benzyl deprotection of phenol")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found and deprotection_found
