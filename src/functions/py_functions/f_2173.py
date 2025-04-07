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
    Detects if the synthesis includes sequential nitrogen functionalization:
    primary amine → secondary amine → tertiary amine
    """
    # Track nitrogen transformations
    primary_to_secondary = False
    secondary_to_tertiary = False

    def dfs_traverse(node):
        nonlocal primary_to_secondary, secondary_to_tertiary

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol:
                    # Check for primary amine to secondary amine transformation
                    primary_amine = Chem.MolFromSmarts("[#7H2]")
                    secondary_amine = Chem.MolFromSmarts("[#7H](-[#6])")
                    tertiary_amine = Chem.MolFromSmarts("[#7](-[#6])(-[#6])-[#6]")

                    # Primary to secondary transformation
                    if any(
                        r and r.HasSubstructMatch(primary_amine) for r in reactants if r
                    ) and product_mol.HasSubstructMatch(secondary_amine):
                        print(
                            "Detected primary amine to secondary amine transformation"
                        )
                        primary_to_secondary = True

                    # Secondary to tertiary transformation
                    if any(
                        r and r.HasSubstructMatch(secondary_amine)
                        for r in reactants
                        if r
                    ) and product_mol.HasSubstructMatch(tertiary_amine):
                        print(
                            "Detected secondary amine to tertiary amine transformation"
                        )
                        secondary_to_tertiary = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if both transformations are detected
    return primary_to_secondary and secondary_to_tertiary
