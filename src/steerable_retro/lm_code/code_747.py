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
    This function detects if the synthesis route involves Boc protection/deprotection strategy.
    """
    boc_pattern = Chem.MolFromSmarts("[CH3]C([CH3])([CH3])[O]C(=O)[N]")
    boc_used = False

    def dfs_traverse(node):
        nonlocal boc_used

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc deprotection: Boc-protected amine -> free amine
                reactant_has_boc = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                        reactant_has_boc = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                product_has_boc = product_mol and product_mol.HasSubstructMatch(boc_pattern)

                if reactant_has_boc and not product_has_boc:
                    print("Boc deprotection detected")
                    boc_used = True
                elif not reactant_has_boc and product_has_boc:
                    print("Boc protection detected")
                    boc_used = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_used
