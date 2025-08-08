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
    This function detects a synthetic strategy involving an isocyanide-based
    multi-component reaction.
    """
    # Initialize tracking variables
    has_isocyanide_mcr = False

    def dfs_traverse(node):
        nonlocal has_isocyanide_mcr

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for multi-component reaction with isocyanide
            if len(reactants_smiles) >= 3:
                isocyanide_pattern = Chem.MolFromSmarts("[#6]-[#7+]#[#6-]")
                has_isocyanide = False

                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(isocyanide_pattern):
                            has_isocyanide = True
                            print("Found isocyanide in MCR")
                            break
                    except:
                        continue

                if has_isocyanide:
                    has_isocyanide_mcr = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_isocyanide_mcr
