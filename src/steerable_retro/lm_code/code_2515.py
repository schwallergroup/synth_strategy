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
    This function detects if the synthesis employs multiple different protecting groups
    (e.g., TBDMS for phenol and THP for heterocycles).
    """
    protecting_groups = {"TBDMS": False, "THP": False}

    def dfs_traverse(node, depth=0):
        nonlocal protecting_groups

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    prod_mol = Chem.MolFromSmiles(product)

                    # Check for TBDMS protection
                    tbdms_patt = Chem.MolFromSmarts("[Si]([C])([C])[C]")
                    if prod_mol.HasSubstructMatch(tbdms_patt):
                        print(f"TBDMS protecting group detected at depth {depth}")
                        protecting_groups["TBDMS"] = True

                    # Check for THP protection
                    thp_patt = Chem.MolFromSmarts("[CH]1[CH2][CH2][CH2][CH2]O1")
                    if prod_mol.HasSubstructMatch(thp_patt):
                        print(f"THP protecting group detected at depth {depth}")
                        protecting_groups["THP"] = True
                except Exception as e:
                    print(f"Error in SMILES processing: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if multiple protecting groups were used
    multiple_pg = sum(protecting_groups.values()) >= 2
    if multiple_pg:
        print(
            f"Multiple protecting groups detected: {[pg for pg, used in protecting_groups.items() if used]}"
        )

    return multiple_pg
