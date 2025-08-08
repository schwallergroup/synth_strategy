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
    This function detects the use of Boc protection/deprotection strategy
    in the synthesis route.
    """
    boc_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal boc_deprotection_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Boc group pattern
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")

                try:
                    # Check if any reactant has Boc group
                    boc_present_in_reactant = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                            boc_present_in_reactant = True
                            break

                    # Check if product doesn't have Boc group
                    product_mol = Chem.MolFromSmiles(product)
                    boc_absent_in_product = product_mol and not product_mol.HasSubstructMatch(
                        boc_pattern
                    )

                    if boc_present_in_reactant and boc_absent_in_product:
                        print("Detected Boc deprotection")
                        boc_deprotection_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_deprotection_detected
