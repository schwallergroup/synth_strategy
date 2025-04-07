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
    This function detects the transformation of a halogen (specifically bromine)
    to a nitrile (CN) group.
    """
    found_br_to_cn = False

    def dfs_traverse(node):
        nonlocal found_br_to_cn

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Br -> CN transformation
            for reactant in reactants:
                if "[Br]" in reactant or "Br" in reactant:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        br_pattern = Chem.MolFromSmarts("[#6][Br]")
                        cn_pattern = Chem.MolFromSmarts("[#6][C]#[N]")

                        if reactant_mol.HasSubstructMatch(
                            br_pattern
                        ) and product_mol.HasSubstructMatch(cn_pattern):
                            # Check if bromine count decreases and nitrile count increases
                            br_count_reactant = len(
                                reactant_mol.GetSubstructMatches(br_pattern)
                            )
                            br_count_product = len(
                                product_mol.GetSubstructMatches(br_pattern)
                            )

                            cn_count_reactant = len(
                                reactant_mol.GetSubstructMatches(cn_pattern)
                            )
                            cn_count_product = len(
                                product_mol.GetSubstructMatches(cn_pattern)
                            )

                            if (
                                br_count_product < br_count_reactant
                                and cn_count_product > cn_count_reactant
                            ):
                                found_br_to_cn = True
                                print("Found bromine to nitrile transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_br_to_cn
