#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthetic route follows a linear synthesis pattern
    (as opposed to convergent synthesis with multiple complex fragments).
    """
    # Track the maximum number of complex reactants in any step
    max_complex_reactants = 0

    def dfs_traverse(node):
        nonlocal max_complex_reactants

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Count complex reactants (those with more than 10 atoms)
                complex_reactant_count = 0
                for r in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r)
                    if r_mol and r_mol.GetNumAtoms() > 10:
                        complex_reactant_count += 1

                max_complex_reactants = max(
                    max_complex_reactants, complex_reactant_count
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If no step has more than one complex reactant, it's a linear synthesis
    is_linear = max_complex_reactants <= 1
    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
