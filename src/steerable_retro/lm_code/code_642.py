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
    This function detects if the synthesis route involves the construction of an imidazopyridine scaffold.
    """
    imidazopyridine_constructed = False
    final_product_has_scaffold = False

    def dfs_traverse(node, depth=0):
        nonlocal imidazopyridine_constructed, final_product_has_scaffold

        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule contains imidazopyridine scaffold
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol is not None:
                    # SMARTS pattern for imidazopyridine scaffold
                    imidazopyridine_pattern = Chem.MolFromSmarts("[n]1[cH][n]c2c1[cH][cH]cc2")
                    if mol.HasSubstructMatch(imidazopyridine_pattern):
                        if depth == 0:  # Final product
                            final_product_has_scaffold = True
                        else:
                            # If an intermediate has the scaffold, it means it was constructed earlier
                            imidazopyridine_constructed = True
            except:
                print("Error processing molecule SMILES for imidazopyridine detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If final product has the scaffold but it wasn't found in intermediates,
    # it means the scaffold was constructed during the synthesis
    if final_product_has_scaffold and not imidazopyridine_constructed:
        imidazopyridine_constructed = True
        print("Imidazopyridine scaffold construction detected")

    return imidazopyridine_constructed
