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
    Detects if the synthesis route maintains a methoxy-substituted aromatic ring throughout.
    """
    all_reactions_have_methoxy_aromatic = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal all_reactions_have_methoxy_aromatic, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for methoxy aromatic in product
            methoxy_aromatic_pattern = Chem.MolFromSmarts("c[O][C]")

            try:
                product_mol = Chem.MolFromSmiles(product)
                if not (product_mol and product_mol.HasSubstructMatch(methoxy_aromatic_pattern)):
                    all_reactions_have_methoxy_aromatic = False
                    print("Found a reaction without methoxy aromatic")
            except:
                all_reactions_have_methoxy_aromatic = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return all_reactions_have_methoxy_aromatic and reaction_count > 0
