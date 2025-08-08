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
    Detects if the synthesis involves thiourea deprotection.
    """
    thiourea_deprotection_detected = False

    def dfs_traverse(node):
        nonlocal thiourea_deprotection_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protected thiourea to free thiourea transformation
                protected_thiourea_pattern = Chem.MolFromSmarts(
                    "[C](=[O])[c]1[cH][cH][cH][cH][cH]1.[NH][C](=[S])[NH]"
                )
                free_thiourea_pattern = Chem.MolFromSmarts("[NH2][C](=[S])[NH]")

                has_protected_thiourea = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(protected_thiourea_pattern):
                        has_protected_thiourea = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    has_protected_thiourea
                    and product_mol
                    and product_mol.HasSubstructMatch(free_thiourea_pattern)
                ):
                    print("Detected thiourea deprotection")
                    thiourea_deprotection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiourea_deprotection_detected
