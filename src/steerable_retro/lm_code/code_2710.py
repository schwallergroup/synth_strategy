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
    This function detects a synthetic strategy involving amine protection/deprotection,
    specifically focusing on Cbz (carboxybenzyl) deprotection.
    """
    has_cbz_deprotection = False

    def dfs_traverse(node):
        nonlocal has_cbz_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Cbz deprotection
            for reactant in reactants:
                react_mol = Chem.MolFromSmiles(reactant)
                if react_mol:
                    # Cbz protected amine pattern
                    cbz_pattern = Chem.MolFromSmarts("[O]=[C]([O][C][c])[N]")
                    if react_mol.HasSubstructMatch(cbz_pattern):
                        # Check if product has primary amine
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            amine_pattern = Chem.MolFromSmarts("[NH2]")
                            if prod_mol.HasSubstructMatch(amine_pattern):
                                has_cbz_deprotection = True
                                print(
                                    f"Detected Cbz deprotection at depth {node['metadata'].get('depth', 'unknown')}"
                                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_cbz_deprotection
