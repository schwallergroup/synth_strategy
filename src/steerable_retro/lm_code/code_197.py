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
    Detects a linear synthesis strategy that utilizes an isothiocyanate intermediate.
    """
    found_isothiocyanate = False
    is_linear_synthesis = True
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal found_isothiocyanate, is_linear_synthesis, reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            reaction_count += 1

            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            max_reactants_per_step = max(max_reactants_per_step, len(reactants_smiles))

            # Check for isothiocyanate
            for smi in reactants_smiles:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    isothiocyanate_pattern = Chem.MolFromSmarts("[#6]-[#7]=[#6]=[#16]")
                    if mol.HasSubstructMatch(isothiocyanate_pattern):
                        found_isothiocyanate = True
                        print("Found isothiocyanate reactant")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A synthesis is considered linear if there are multiple reactions and
    # no more than 2 reactants per step on average
    is_linear_synthesis = reaction_count > 1 and max_reactants_per_step <= 2

    return found_isothiocyanate and is_linear_synthesis
