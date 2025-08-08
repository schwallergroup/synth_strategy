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
    Detects if the synthesis incorporates an N-methylpiperazine fragment
    in a late-stage nucleophilic substitution.
    """
    has_n_methylpiperazine = False
    has_late_stage_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_methylpiperazine, has_late_stage_incorporation

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check for N-methylpiperazine
            n_methylpiperazine_pattern = Chem.MolFromSmarts("[CH3][N]1[CH2][CH2][N][CH2][CH2]1")
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(n_methylpiperazine_pattern):
                        has_n_methylpiperazine = True
                        # If depth is 0 or 1, it's considered late-stage
                        if depth <= 1:
                            has_late_stage_incorporation = True
                            print(f"Found N-methylpiperazine at depth {depth} (late-stage)")
                        else:
                            print(f"Found N-methylpiperazine at depth {depth}")
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if N-methylpiperazine is incorporated in a late stage
    return has_n_methylpiperazine and has_late_stage_incorporation
