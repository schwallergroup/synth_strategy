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
    Detects a synthetic sequence involving nitro reduction to amine
    followed by tetrazole formation from the amine.
    """
    # Track if we've found the key transformations
    found_nitro_reduction = False
    found_tetrazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction, found_tetrazole_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction
            if not found_nitro_reduction:
                nitro_patt = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_patt = Chem.MolFromSmarts("c-[NH2]")

                reactant_has_nitro = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(nitro_patt):
                        reactant_has_nitro = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_nitro
                    and product_mol
                    and product_mol.HasSubstructMatch(amine_patt)
                ):
                    print("Found nitro reduction to amine")
                    found_nitro_reduction = True

            # Check for tetrazole formation from amine
            if not found_tetrazole_formation:
                amine_patt = Chem.MolFromSmarts("c-[NH2]")
                tetrazole_patt = Chem.MolFromSmarts("c-n1nnnc1C")

                reactant_has_amine = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(amine_patt):
                        reactant_has_amine = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                if (
                    reactant_has_amine
                    and product_mol
                    and product_mol.HasSubstructMatch(tetrazole_patt)
                ):
                    print("Found tetrazole formation from amine")
                    found_tetrazole_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if both transformations were found
    return found_nitro_reduction and found_tetrazole_formation
