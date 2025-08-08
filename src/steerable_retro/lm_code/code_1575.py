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
    Detects if certain functional groups (carboxylic acid, hydroxyl) are preserved
    throughout the synthesis.
    """
    # Track functional groups at each depth
    depth_to_functional_groups = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and not node.get("in_stock", False):
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for carboxylic acid
                carboxylic_acid = mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]"))
                # Check for hydroxyl (not part of carboxylic acid)
                hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))

                if depth not in depth_to_functional_groups:
                    depth_to_functional_groups[depth] = {
                        "carboxylic_acid": carboxylic_acid,
                        "hydroxyl": hydroxyl,
                    }
                else:
                    depth_to_functional_groups[depth]["carboxylic_acid"] |= carboxylic_acid
                    depth_to_functional_groups[depth]["hydroxyl"] |= hydroxyl

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if functional groups are preserved across at least 3 depths
    if len(depth_to_functional_groups) >= 3:
        carboxylic_preserved = all(
            info["carboxylic_acid"] for info in depth_to_functional_groups.values()
        )
        hydroxyl_preserved = all(info["hydroxyl"] for info in depth_to_functional_groups.values())

        if carboxylic_preserved:
            print("Carboxylic acid group preserved throughout synthesis")
        if hydroxyl_preserved:
            print("Hydroxyl group preserved throughout synthesis")

        return carboxylic_preserved or hydroxyl_preserved

    return False
