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
    Detects if the synthesis involves protection of an amine with a carbamate group
    (specifically benzyl carbamate/Cbz) that is carried through multiple steps.
    """
    carbamate_formation = False
    protected_steps = 0

    def dfs_traverse(node):
        nonlocal carbamate_formation, protected_steps

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carbamate formation
                amine_pattern = Chem.MolFromSmarts("[#7;H2]")
                carbamate_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")

                # Check if this reaction forms a carbamate from an amine
                has_amine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                        break

                prod_mol = Chem.MolFromSmiles(product)
                if has_amine and prod_mol and prod_mol.HasSubstructMatch(carbamate_pattern):
                    carbamate_formation = True
                    print("Found carbamate protection of amine")

                # Count steps where protected amine is carried through
                if prod_mol and prod_mol.HasSubstructMatch(carbamate_pattern):
                    protected_steps += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if carbamate is formed and carried through at least 2 more steps
    return carbamate_formation and protected_steps >= 3
