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
    This function detects if the synthesis includes an amide bond formation step.
    """
    amide_pattern = Chem.MolFromSmarts("[#6]C(=O)[#7]")
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain carboxylic acid or derivative and amine
            acid_pattern = Chem.MolFromSmarts("[#6]C(=O)[#8]")
            amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

            product_mol = Chem.MolFromSmiles(product)

            has_acid = False
            has_amine = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    if reactant_mol.HasSubstructMatch(acid_pattern):
                        has_acid = True
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            if (
                has_acid
                and has_amine
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                amide_formation_found = True
                print("Amide bond formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
