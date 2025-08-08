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
    This function detects a synthetic strategy involving conversion of a dibromo-alkene
    to a terminal bromoacetylene.
    """
    dibromo_to_alkyne_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal dibromo_to_alkyne_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains terminal bromoacetylene
                product_mol = Chem.MolFromSmiles(product)
                bromo_alkyne_pattern = Chem.MolFromSmarts("[#6]#[#6][Br]")

                if product_mol and product_mol.HasSubstructMatch(bromo_alkyne_pattern):
                    # Check if reactant contains dibromo-alkene
                    dibromo_alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]([Br])[Br]")
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(dibromo_alkene_pattern):
                            print(
                                f"Dibromo-alkene to bromoacetylene conversion detected at depth {depth}"
                            )
                            dibromo_to_alkyne_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return dibromo_to_alkyne_detected
