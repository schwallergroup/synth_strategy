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
    This function detects if a nitrile group is carried through the synthesis without modification.
    """
    nitrile_reactions = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for nitrile in product
                nitrile_pattern = Chem.MolFromSmarts("[#6]-C#N")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_reactions.append(node.get("metadata", {}).get("ID", "unknown"))

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Found nitrile group in {len(nitrile_reactions)} reactions")
    return len(nitrile_reactions) >= 1  # At least one reaction involves a nitrile
