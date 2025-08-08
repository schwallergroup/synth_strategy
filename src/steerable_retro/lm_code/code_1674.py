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
    Detects if the synthesis follows a linear strategy (no convergent steps with multiple complex fragments).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                # In retrosynthetic analysis, we look at the products of the forward reaction
                # which are the reactants in our retrosynthetic direction
                reactants = rsmi.split(">")[2].split(".")

                # Count complex reactants (more than 10 heavy atoms)
                complex_reactants = 0
                print(f"Examining reaction: {rsmi}")
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        heavy_atoms = mol.GetNumHeavyAtoms()
                        print(f"  Reactant: {reactant}, Heavy atoms: {heavy_atoms}")
                        if heavy_atoms > 10:
                            complex_reactants += 1

                # If more than one complex reactant and multiple reactants, it's a convergent step
                if complex_reactants > 1 and len(reactants) > 1:
                    is_linear = False
                    print(f"Convergent step detected with {complex_reactants} complex reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
