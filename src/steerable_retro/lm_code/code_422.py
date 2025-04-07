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
    This function detects if the synthesis follows a convergent approach by identifying
    a reaction that combines two complex fragments.
    """
    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if we have at least two complex reactants
            complex_reactants = 0
            indole_found = False
            fluorinated_ring_found = False

            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                # Check for indole core
                indole_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#7][#6]2[#6]([#6]1)[#6][#6][#6][#6]2"
                )
                if reactant_mol.HasSubstructMatch(indole_pattern):
                    indole_found = True
                    complex_reactants += 1

                # Check for fluorinated aromatic ring
                fluoro_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]([F])[#6][#6][#6]1")
                if reactant_mol.HasSubstructMatch(fluoro_pattern):
                    fluorinated_ring_found = True
                    complex_reactants += 1

                # Count atoms as a measure of complexity
                if reactant_mol.GetNumAtoms() > 10:
                    complex_reactants += 1

            # If we have at least two complex reactants or both key fragments
            if complex_reactants >= 2 or (indole_found and fluorinated_ring_found):
                is_convergent = True
                print(f"Convergent synthesis detected at depth {depth}")
                print(
                    f"Found indole: {indole_found}, Found fluorinated ring: {fluorinated_ring_found}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Convergent synthesis strategy: {is_convergent}")
    return is_convergent
