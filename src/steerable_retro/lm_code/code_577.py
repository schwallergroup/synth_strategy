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
    Detects if the synthesis involves an ether formation step.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol/phenol in reactants
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                # Check for ether in product
                ether_pattern = Chem.MolFromSmarts("[#6]-[O]-[#6]")

                has_alcohol = sum(
                    1
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                )
                has_ether = (
                    Chem.MolFromSmiles(product).HasSubstructMatch(ether_pattern)
                    if Chem.MolFromSmiles(product)
                    else False
                )

                if has_alcohol >= 1 and has_ether:
                    print("Detected ether formation step")
                    result = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
