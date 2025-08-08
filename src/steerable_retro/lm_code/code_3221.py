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
    Detects if the synthesis involves a thiol deprotection step.
    """
    thiol_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal thiol_deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiol in product
                thiol_pattern = Chem.MolFromSmarts("[#6]-[#16H]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(thiol_pattern):
                    # Check if reactants don't have free thiol
                    has_free_thiol_in_reactants = False

                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(thiol_pattern):
                            has_free_thiol_in_reactants = True

                    if not has_free_thiol_in_reactants:
                        print("Thiol deprotection detected")
                        thiol_deprotection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiol_deprotection_detected
