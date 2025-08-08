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
    Detects if the synthesis route uses a convergent approach with multiple fragments.
    Looks for reactions with multiple complex reactants (3+ heavy atoms each).
    """
    has_convergent = False

    def dfs_traverse(node):
        nonlocal has_convergent

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count complex reactants (3+ heavy atoms)
            complex_reactants = 0
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                heavy_atom_count = reactant_mol.GetNumHeavyAtoms()
                if heavy_atom_count >= 3:
                    complex_reactants += 1

            # If 2+ complex reactants, it's a convergent step
            if complex_reactants >= 2:
                has_convergent = True
                print(f"Found convergent synthesis step with {complex_reactants} complex reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_convergent
