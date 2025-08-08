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
    Detects a strategy involving a trifluoromethyl ether group maintained throughout the synthesis.
    """
    # Track if we found the trifluoromethyl ether in all molecules
    all_molecules_have_tfm_ether = True
    molecule_count = 0

    def dfs_traverse(node):
        nonlocal all_molecules_have_tfm_ether, molecule_count

        if node["type"] == "mol" and not node.get("in_stock", False):
            molecule_count += 1
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for trifluoromethyl ether pattern: O-CH2-CF3
                tfm_ether_patt = Chem.MolFromSmarts("[O][C][C]([F])([F])[F]")
                if not mol.HasSubstructMatch(tfm_ether_patt):
                    all_molecules_have_tfm_ether = False
                    print(f"Molecule without trifluoromethyl ether: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if all non-starting molecules have the trifluoromethyl ether group
    # and we've examined at least one molecule
    return all_molecules_have_tfm_ether and molecule_count > 0
