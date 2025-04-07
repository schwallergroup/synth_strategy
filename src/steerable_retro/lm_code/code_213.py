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
    This function detects if the synthesis route involves amide formation
    in the late stages of the synthesis (low depth).
    """
    late_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[O])-[#8;H1]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")
                amide_pattern = Chem.MolFromSmarts("[#6](=[O])-[#7]")

                has_acid = False
                has_amine = False

                for r in reactants:
                    if r:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(carboxylic_acid_pattern):
                                has_acid = True
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                # Check for amide in product
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and has_acid
                    and has_amine
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    print(f"Late-stage amide formation detected at depth {depth}")
                    late_amide_formation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_amide_formation
