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
    Detects if the synthesis follows a linear strategy rather than a convergent one.
    Linear synthesis typically has only one complex building block per reaction.
    """
    linear_strategy = True
    complex_reactant_threshold = (
        15  # Minimum atom count to consider a reactant "complex"
    )

    def dfs_traverse(node):
        nonlocal linear_strategy

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (those with more than threshold atoms)
                complex_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > complex_reactant_threshold:
                            complex_reactants += 1
                    except:
                        continue

                # If more than one complex reactant, it's likely a convergent step
                if complex_reactants > 1:
                    print(
                        f"Detected convergent step with {complex_reactants} complex reactants"
                    )
                    linear_strategy = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Synthesis strategy is {'linear' if linear_strategy else 'convergent'}")
    return linear_strategy
