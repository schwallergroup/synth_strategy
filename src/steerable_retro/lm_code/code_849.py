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
    This function detects if the synthetic route contains amide bond formation reactions.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride pattern
            acid_chloride_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[Cl]")

            # Check for amine pattern
            amine_pattern = Chem.MolFromSmarts("[#6][NH2]")

            # Check for amide pattern in product
            amide_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[#7][#6]")

            # Check if reactants contain acid chloride and amine, and product contains amide
            acid_chloride_present = False
            amine_present = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                    acid_chloride_present = True

                if reactant_mol.HasSubstructMatch(amine_pattern):
                    amine_present = True

            product_mol = Chem.MolFromSmiles(product)
            if (
                acid_chloride_present
                and amine_present
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                amide_formation_found = True
                print("Found amide bond formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return amide_formation_found
