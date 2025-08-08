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
    This function detects a synthetic strategy involving multiple heteroatom
    bond formations (C-N, C-O) throughout the synthesis.
    """
    c_n_bond_formations = 0
    c_o_bond_formations = 0

    def dfs_traverse(node):
        nonlocal c_n_bond_formations, c_o_bond_formations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Create molecules for analysis
            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check for C-N bond formation
                c_n_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                c_n_in_reactants = len(reactants_mol.GetSubstructMatches(c_n_pattern))
                c_n_in_product = len(product_mol.GetSubstructMatches(c_n_pattern))

                if c_n_in_product > c_n_in_reactants:
                    c_n_bond_formations += 1
                    print(
                        f"C-N bond formation detected: {c_n_in_product - c_n_in_reactants} new bonds"
                    )

                # Check for C-O bond formation
                c_o_pattern = Chem.MolFromSmarts("[#6]-[#8]")
                c_o_in_reactants = len(reactants_mol.GetSubstructMatches(c_o_pattern))
                c_o_in_product = len(product_mol.GetSubstructMatches(c_o_pattern))

                if c_o_in_product > c_o_in_reactants:
                    c_o_bond_formations += 1
                    print(
                        f"C-O bond formation detected: {c_o_in_product - c_o_in_reactants} new bonds"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if at least one C-N and one C-O bond formation is detected
    return c_n_bond_formations > 0 and c_o_bond_formations > 0
