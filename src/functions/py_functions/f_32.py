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
    This function detects if the synthesis involves trifluoroacetyl protection
    of an amine.
    """
    trifluoroacetyl_protection_detected = False

    def dfs_traverse(node):
        nonlocal trifluoroacetyl_protection_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for trifluoroacetyl protection: NH2 + trifluoroacetic derivative -> NH-COCF3
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            trifluoroacetyl_pattern = Chem.MolFromSmarts(
                "[#7]-[#6](=[#8])-[#6]([#9])([#9])[#9]"
            )

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            has_amine = any(
                r and r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
            )

            if (
                has_amine
                and product_mol
                and product_mol.HasSubstructMatch(trifluoroacetyl_pattern)
            ):
                print("Detected trifluoroacetyl protection of amine")
                trifluoroacetyl_protection_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return trifluoroacetyl_protection_detected
