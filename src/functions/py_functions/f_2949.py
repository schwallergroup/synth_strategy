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
    This function detects if the synthetic route employs a convergent approach
    where two complex fragments are joined.
    """
    convergent_step_found = False

    def dfs_traverse(node):
        nonlocal convergent_step_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a coupling reaction with at least 2 reactants
            if len(reactants) >= 2:
                # Check complexity of reactants (number of atoms as a simple metric)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                reactant_complexities = [
                    mol.GetNumAtoms() if mol is not None else 0 for mol in reactant_mols
                ]

                # Consider it convergent if both main reactants have significant complexity
                # (arbitrary threshold of 10 atoms)
                complex_reactants = sum(
                    1 for complexity in reactant_complexities if complexity >= 10
                )

                if complex_reactants >= 2:
                    print(
                        f"Found convergent synthesis step with {complex_reactants} complex reactants"
                    )
                    convergent_step_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return convergent_step_found
