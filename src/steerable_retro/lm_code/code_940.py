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
    This function detects the use of TBDMS (tert-butyldimethylsilyl) protection
    of alcohols in the synthetic route.
    """
    tbdms_pattern = Chem.MolFromSmarts("[#6]-[#14](-[#6])(-[#6])-[#8]-[#6]")
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains TBDMS group but reactants don't
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(tbdms_pattern):
                    # Check if this is a protection reaction (OH → OTBDMS)
                    reactant_has_oh = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#8H1]-[#6]")
                        ):
                            reactant_has_oh = True
                            break

                    if reactant_has_oh:
                        protection_found = True
                        print("TBDMS protection detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protection_found
