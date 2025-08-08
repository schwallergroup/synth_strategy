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
    This function detects if the synthesis involves sequential formation of different heteroatom bonds.
    """
    # Track bond formations at different depths
    bond_formations = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)

                # Check for different bond formations
                c_n_amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                c_o_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")
                s_n_sulfonamide_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[N]")
                c_s_pattern = Chem.MolFromSmarts("[c][S]")

                if product_mol:
                    if product_mol.HasSubstructMatch(c_n_amide_pattern):
                        bond_formations[depth] = "C-N amide"
                    elif product_mol.HasSubstructMatch(c_o_ether_pattern):
                        bond_formations[depth] = "C-O ether"
                    elif product_mol.HasSubstructMatch(s_n_sulfonamide_pattern):
                        bond_formations[depth] = "S-N sulfonamide"
                    elif product_mol.HasSubstructMatch(c_s_pattern):
                        bond_formations[depth] = "C-S"

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if we have at least 3 different types of heteroatom bond formations
    unique_bond_types = set(bond_formations.values())
    print(f"Detected bond formations: {bond_formations}")
    print(f"Unique bond types: {unique_bond_types}")

    return len(unique_bond_types) >= 3
