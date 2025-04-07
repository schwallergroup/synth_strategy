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
    Detects if the synthesis route employs an early-stage benzylation strategy,
    where a benzyl group is attached to a thiophene core.
    """
    benzylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal benzylation_detected

        if node["type"] == "reaction" and depth >= 3:  # Focus on early-stage reactions
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiophene in reactants
                thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#16][#6]1")
                # Check for benzyl bromide in reactants
                benzyl_bromide_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1[C][Br]")
                # Check for benzylated thiophene in product
                benzylated_thiophene_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#16][#6]1[C][#6]2[#6][#6][#6][#6][#6]2"
                )

                thiophene_present = False
                benzyl_bromide_present = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(thiophene_pattern):
                            thiophene_present = True
                        if mol.HasSubstructMatch(benzyl_bromide_pattern):
                            benzyl_bromide_present = True

                product_mol = Chem.MolFromSmiles(product)
                if (
                    thiophene_present
                    and benzyl_bromide_present
                    and product_mol
                    and product_mol.HasSubstructMatch(benzylated_thiophene_pattern)
                ):
                    print(f"Early-stage benzylation detected at depth {depth}")
                    benzylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return benzylation_detected
