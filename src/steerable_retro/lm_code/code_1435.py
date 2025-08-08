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
    This function detects if the synthetic route involves nitrogen-rich compounds
    with multiple nitrogen-containing functional groups.
    """
    nitrogen_count = 0
    n_functional_groups = 0

    def dfs_traverse(node):
        nonlocal nitrogen_count, n_functional_groups

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count nitrogen atoms
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        nitrogen_count += 1

                # Check for nitrogen-containing functional groups
                n_functional_groups_patterns = [
                    Chem.MolFromSmarts("[NX3]"),  # amine
                    Chem.MolFromSmarts("[NX3][NX3]"),  # hydrazine
                    Chem.MolFromSmarts("[#7]=[#6]=[#8]"),  # isocyanate
                    Chem.MolFromSmarts("[NX3][CX3]=[O]"),  # amide
                    Chem.MolFromSmarts("[NX3+](=[O])[O-]"),  # nitro
                    Chem.MolFromSmarts("[n]"),  # aromatic nitrogen
                ]

                for pattern in n_functional_groups_patterns:
                    if mol.HasSubstructMatch(pattern):
                        n_functional_groups += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Consider it nitrogen-rich if it has at least 3 nitrogen atoms and 2 different N-functional groups
    is_nitrogen_rich = nitrogen_count >= 3 and n_functional_groups >= 2
    if is_nitrogen_rich:
        print(
            f"Nitrogen-rich synthesis detected: {nitrogen_count} N atoms, {n_functional_groups} N-functional groups"
        )

    return is_nitrogen_rich
