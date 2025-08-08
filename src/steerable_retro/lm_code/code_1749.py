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
    This function detects Boc deprotection to form a primary amine in the early stages of synthesis.
    """
    boc_deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_deprotection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc group in reactants
            boc_pattern = Chem.MolFromSmarts("C(=O)OC(C)(C)C")
            nh2_pattern = Chem.MolFromSmarts("[NH2][#6]")

            for reactant in reactants:
                try:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol and r_mol.HasSubstructMatch(boc_pattern):
                        # Check if product has primary amine
                        p_mol = Chem.MolFromSmiles(product)
                        if p_mol and p_mol.HasSubstructMatch(nh2_pattern):
                            print("Found Boc deprotection to primary amine")
                            boc_deprotection_found = True
                except:
                    continue

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_deprotection_found
