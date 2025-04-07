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
    Detects if the synthesis route involves a sequential transformation
    from aldehyde to oxime to amine.
    """
    aldehyde_depths = []
    oxime_depths = []
    amine_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal aldehyde_depths, oxime_depths, amine_depths

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
            oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2H]")
            amine_pattern = Chem.MolFromSmarts("[CX4][NX3]")

            product_mol = Chem.MolFromSmiles(product)

            # Check for aldehyde to oxime transformation
            reactants_have_aldehyde = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(aldehyde_pattern):
                    reactants_have_aldehyde = True
                    aldehyde_depths.append(depth)
                    print(f"Aldehyde detected in reactant at depth {depth}")
                    break

            if (
                reactants_have_aldehyde
                and product_mol
                and product_mol.HasSubstructMatch(oxime_pattern)
            ):
                oxime_depths.append(depth)
                print(f"Oxime formation from aldehyde detected at depth {depth}")

            # Check for oxime to amine transformation
            reactants_have_oxime = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(oxime_pattern):
                    reactants_have_oxime = True
                    if depth not in oxime_depths:
                        oxime_depths.append(depth)
                        print(f"Oxime detected in reactant at depth {depth}")
                    break

            if (
                reactants_have_oxime
                and product_mol
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                amine_depths.append(depth)
                print(f"Amine formation from oxime detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have the complete sequence
    if aldehyde_depths and oxime_depths and amine_depths:
        # Check if there's a sequential relationship
        for aldehyde_depth in aldehyde_depths:
            for oxime_depth in oxime_depths:
                for amine_depth in amine_depths:
                    if (
                        amine_depth <= oxime_depth <= aldehyde_depth
                    ):  # Remember lower depth means later in synthesis
                        print(f"Complete aldehyde→oxime→amine sequence detected")
                        return True

    return False
