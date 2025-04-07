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
    This function detects if the synthesis uses late-stage amide formation
    as the final or near-final step.
    """
    found_late_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amide

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        # Check if reactants contain acid and amine
                        has_acid = False
                        has_amine = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("[C](=[O])[O]")
                                ):
                                    has_acid = True
                                if reactant_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("[N]")
                                ):
                                    has_amine = True

                        if has_acid and has_amine:
                            print(f"Found late-stage amide formation at depth {depth}")
                            found_late_amide = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_late_amide
