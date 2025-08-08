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
    Detects if the synthesis involves compounds containing fluorophenyl groups
    that are maintained throughout the synthesis.
    """
    final_product_has_fluorophenyl = False
    intermediates_with_fluorophenyl = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_fluorophenyl, intermediates_with_fluorophenyl

        if node["type"] == "mol":
            # Check for fluorophenyl groups
            fluorophenyl_pattern = Chem.MolFromSmarts("c1cc(F)ccc1")

            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(fluorophenyl_pattern):
                    if depth == 0:  # Final product
                        final_product_has_fluorophenyl = True
                        print("Final product contains fluorophenyl group")
                    else:
                        intermediates_with_fluorophenyl += 1
                        print(f"Intermediate at depth {depth} contains fluorophenyl group")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return final_product_has_fluorophenyl and intermediates_with_fluorophenyl > 0
