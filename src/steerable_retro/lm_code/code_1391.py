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
    This function detects a protection/deprotection sequence in the synthesis,
    specifically looking for Cbz protection of amines.
    """
    has_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_deprotection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Cbz deprotection
            if len(reactants) == 1:
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                # Pattern for Cbz-protected amine
                cbz_pattern = Chem.MolFromSmarts("[O]=[C][O][C][c]1[cH][cH][cH][cH][cH]1")
                # Pattern for amine
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                if (
                    reactant_mol.HasSubstructMatch(cbz_pattern)
                    and product_mol.HasSubstructMatch(amine_pattern)
                    and not product_mol.HasSubstructMatch(cbz_pattern)
                ):
                    has_deprotection = True
                    print("Found Cbz deprotection")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_deprotection
