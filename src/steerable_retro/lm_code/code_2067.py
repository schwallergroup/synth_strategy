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
    This function detects if the synthesis follows a linear strategy
    rather than a convergent one.
    """
    # In a linear synthesis, most reactions have only one non-reagent reactant
    # We'll count reactions and check how many have multiple complex reactants

    total_reactions = 0
    convergent_reactions = 0
    threshold_complex_atoms = 10  # Define complex molecule as having > 10 atoms

    def dfs_traverse(node):
        nonlocal total_reactions, convergent_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants
                complex_reactants = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > threshold_complex_atoms:
                        complex_reactants += 1

                if complex_reactants > 1:
                    convergent_reactions += 1
                    print(f"Convergent reaction detected: {complex_reactants} complex reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If less than 25% of reactions are convergent, consider it a linear synthesis
    is_linear = total_reactions > 0 and (convergent_reactions / total_reactions) < 0.25
    if is_linear:
        print(
            f"Linear synthesis detected: {convergent_reactions}/{total_reactions} convergent reactions"
        )

    return is_linear
