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
    Detects amide formation from acid chloride and amine.
    """
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for acid chloride and amine patterns in reactants
            has_acid_chloride = False
            has_amine = False

            for reactant in reactants:
                r_mol = Chem.MolFromSmiles(reactant)
                if not r_mol:
                    continue

                acid_chloride_patt = Chem.MolFromSmarts("[C](=O)[Cl]")
                if r_mol.HasSubstructMatch(acid_chloride_patt):
                    has_acid_chloride = True

                amine_patt = Chem.MolFromSmarts("[NH2]")
                if r_mol.HasSubstructMatch(amine_patt):
                    has_amine = True

            # Check for amide pattern in product
            p_mol = Chem.MolFromSmiles(product)
            has_amide = False
            if p_mol:
                amide_patt = Chem.MolFromSmarts("[C](=O)[NH]")
                if p_mol.HasSubstructMatch(amide_patt):
                    has_amide = True

            if has_acid_chloride and has_amine and has_amide:
                has_amide_formation = True
                print("Detected amide formation from acid chloride and amine")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_amide_formation
