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
    Detects a late-stage O-alkylation where a phenol is alkylated with a hydroxyethyl chain
    in the final steps of the synthesis.
    """
    o_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_found

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at depth 0 or 1 (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for O-alkylation pattern
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                hydroxyethyl_ether_pattern = Chem.MolFromSmarts("[c][O][C][C][OH]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and product_mol.HasSubstructMatch(
                    hydroxyethyl_ether_pattern
                ):
                    phenol_found = any(
                        mol and mol.HasSubstructMatch(phenol_pattern)
                        for mol in reactant_mols
                    )

                    if phenol_found:
                        o_alkylation_found = True
                        print(f"Found late-stage O-alkylation at depth {depth}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return o_alkylation_found
