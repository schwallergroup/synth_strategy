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
    This function detects if the synthesis route involves incorporation of a thiophene ring.
    """
    thiophene_incorporation = False

    def dfs_traverse(node):
        nonlocal thiophene_incorporation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiophene in reactants
            thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#16][#6][#6]1")
            thiophene_in_reactants = False

            for reactant in reactants:
                try:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(thiophene_pattern):
                        thiophene_in_reactants = True
                        break
                except:
                    pass

            # Check for thiophene in product
            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    thiophene_in_reactants
                    and product_mol
                    and product_mol.HasSubstructMatch(thiophene_pattern)
                ):
                    print("Thiophene incorporation detected")
                    thiophene_incorporation = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return thiophene_incorporation
