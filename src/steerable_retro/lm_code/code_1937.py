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
    This function detects if the synthetic route follows a linear synthesis strategy
    (as opposed to convergent).
    """
    # In a linear synthesis, most reactions have only one non-commercial reactant
    # We'll count reactions and check how many have multiple complex reactants
    total_reactions = 0
    convergent_reactions = 0

    def is_complex_molecule(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Define criteria for "complex" molecules (e.g., more than 10 atoms)
            return mol.GetNumAtoms() > 10
        return False

    def dfs_traverse(node):
        nonlocal total_reactions, convergent_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                total_reactions += 1

                # Count complex reactants
                complex_reactants = sum(1 for r in reactants_smiles if is_complex_molecule(r))

                if complex_reactants > 1:
                    convergent_reactions += 1
                    print(f"Convergent reaction detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If more than 25% of reactions are convergent, we don't consider it a linear strategy
    if total_reactions > 0:
        return (convergent_reactions / total_reactions) <= 0.25

    return False
