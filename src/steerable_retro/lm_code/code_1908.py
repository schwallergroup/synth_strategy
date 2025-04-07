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
    This function detects if the synthesis uses a linear assembly strategy
    rather than a convergent approach.
    """
    # In a linear synthesis, most reactions have only one non-reagent reactant
    # We'll count reactions and check how many have multiple complex reactants

    total_reactions = 0
    convergent_reactions = 0

    def is_complex_molecule(smiles):
        # Define a complex molecule as having more than 10 atoms or multiple rings
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if mol.GetNumAtoms() > 10:
                    return True
                ring_info = mol.GetRingInfo()
                if ring_info.NumRings() > 1:
                    return True
        except:
            pass
        return False

    def dfs_traverse(node):
        nonlocal total_reactions, convergent_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            total_reactions += 1

            # Count complex reactants
            complex_reactant_count = sum(1 for r in reactants if is_complex_molecule(r))

            if complex_reactant_count >= 2:
                convergent_reactions += 1
                print(f"Found convergent reaction: {rsmi}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If less than 25% of reactions are convergent, consider it a linear strategy
    is_linear = total_reactions > 0 and (convergent_reactions / total_reactions) < 0.25

    print(f"Total reactions: {total_reactions}, Convergent reactions: {convergent_reactions}")
    print(f"Uses linear fragment assembly strategy: {is_linear}")

    return is_linear
