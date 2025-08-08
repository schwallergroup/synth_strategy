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
    Detects a convergent synthesis involving a piperazine scaffold where
    two complex fragments are combined in a late-stage reaction.
    """
    # Track if we found a convergent step
    convergent_step_found = False
    # Track if piperazine is present
    piperazine_present = False

    def dfs_traverse(node):
        nonlocal convergent_step_found, piperazine_present

        if node["type"] == "mol":
            # Check if molecule contains piperazine
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[N]1[C][C][N][C][C]1")):
                    piperazine_present = True

        elif node["type"] == "reaction":
            # Check if this is a convergent step (combining multiple fragments)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If we have multiple complex reactants, consider it convergent
                complex_reactants = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:  # Arbitrary threshold for "complex"
                        complex_reactants += 1

                if complex_reactants >= 2:
                    print(f"Found convergent step with {complex_reactants} complex reactants")
                    convergent_step_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return convergent_step_found and piperazine_present
