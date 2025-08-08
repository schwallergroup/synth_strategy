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
    Detects if the synthesis route includes phenol alkylation with a dibromoalkane
    """
    phenol_alkylation_found = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if one reactant is a dibromoalkane
            dibromo_pattern = Chem.MolFromSmarts("[Br][C][C][C][Br]")
            phenol_pattern = Chem.MolFromSmarts("[c][OH]")
            aryl_ether_pattern = Chem.MolFromSmarts("[c][O][C][C][C][Br]")

            dibromo_found = False
            phenol_found = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(dibromo_pattern):
                        dibromo_found = True
                    if reactant_mol.HasSubstructMatch(phenol_pattern):
                        phenol_found = True

            product_mol = Chem.MolFromSmiles(product)

            if (
                dibromo_found
                and phenol_found
                and product_mol
                and product_mol.HasSubstructMatch(aryl_ether_pattern)
            ):
                print("Found phenol alkylation with dibromoalkane")
                phenol_alkylation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return phenol_alkylation_found
