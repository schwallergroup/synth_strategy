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
    This function detects convergent synthesis with amide bond formation between two complex fragments.
    """
    convergent_amide_found = False

    def dfs_traverse(node):
        nonlocal convergent_amide_found

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride and amine patterns in reactants
            acid_chloride_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[Cl]")
            amine_pattern = Chem.MolFromSmarts("[NX3;H2]")
            amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;H1]")

            # Check if both reactants are complex (more than 15 atoms)
            complex_reactants = 0
            acid_chloride_found = False
            amine_found = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.GetNumAtoms() > 15:
                        complex_reactants += 1
                    if mol.HasSubstructMatch(acid_chloride_pattern):
                        acid_chloride_found = True
                    if mol.HasSubstructMatch(amine_pattern):
                        amine_found = True

            product_mol = Chem.MolFromSmiles(product)
            amide_formed = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            if (
                acid_chloride_found
                and amine_found
                and amide_formed
                and complex_reactants >= 1
            ):
                convergent_amide_found = True
                print("Found convergent synthesis with amide bond formation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return convergent_amide_found
