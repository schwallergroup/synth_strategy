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
    This function detects a convergent synthesis strategy where two complex fragments
    are joined in the middle of the synthesis.
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if we have multiple complex reactants
            # Complexity can be estimated by molecule size, number of rings, etc.
            if len(reactants) >= 2:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                # Count atoms as a simple complexity measure
                atom_counts = [mol.GetNumAtoms() if mol else 0 for mol in reactant_mols]

                # If we have at least two reactants with significant complexity
                if sum(count >= 10 for count in atom_counts) >= 2:
                    print(f"Complex fragment joining detected at depth {depth}")
                    # Middle of synthesis (not first or last step)
                    if depth > 0:
                        convergent_synthesis_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_synthesis_detected
