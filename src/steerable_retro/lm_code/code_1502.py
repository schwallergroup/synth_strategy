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
    Detects if the synthetic route contains a tertiary amine fragment,
    particularly N,N-dimethylaminoethyl groups.
    """
    has_tertiary_amine = False

    def dfs_traverse(node):
        nonlocal has_tertiary_amine

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Pattern for tertiary amine, specifically N,N-dimethylaminoethyl
                tertiary_amine_pattern = Chem.MolFromSmarts("[#6]-[#7]([#6])[#6]")
                dimethylaminoethyl_pattern = Chem.MolFromSmarts("[#6]-[#6]-[#7]([#6])[#6]")

                if mol.HasSubstructMatch(dimethylaminoethyl_pattern):
                    has_tertiary_amine = True
                    print(f"Found tertiary amine fragment in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_tertiary_amine
