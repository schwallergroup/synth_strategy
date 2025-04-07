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
    Detects a synthetic strategy involving late-stage SNAr coupling to form diaryl ether
    followed by benzofuran ring formation.
    """
    # Track if we found the key transformations
    found_diaryl_ether_formation = False
    found_benzofuran_formation = False
    diaryl_ether_depth = -1
    benzofuran_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_diaryl_ether_formation, found_benzofuran_formation
        nonlocal diaryl_ether_depth, benzofuran_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and all(r is not None for r in reactants):
                    # Check for diaryl ether formation (SNAr)
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                    fluoro_pattern = Chem.MolFromSmarts("[c][F,Cl,Br,I]")
                    diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

                    has_phenol = any(
                        r.HasSubstructMatch(phenol_pattern) for r in reactants
                    )
                    has_haloarene = any(
                        r.HasSubstructMatch(fluoro_pattern) for r in reactants
                    )
                    product_has_diaryl_ether = product.HasSubstructMatch(
                        diaryl_ether_pattern
                    )

                    if has_phenol and has_haloarene and product_has_diaryl_ether:
                        found_diaryl_ether_formation = True
                        diaryl_ether_depth = depth
                        print(f"Found diaryl ether formation at depth {depth}")

                    # Check for benzofuran formation
                    benzofuran_pattern = Chem.MolFromSmarts(
                        "[c]1[c][c][c][c]2[c]1[o][c][c]2"
                    )

                    reactants_have_benzofuran = any(
                        r.HasSubstructMatch(benzofuran_pattern) for r in reactants
                    )
                    product_has_benzofuran = product.HasSubstructMatch(
                        benzofuran_pattern
                    )

                    if not reactants_have_benzofuran and product_has_benzofuran:
                        found_benzofuran_formation = True
                        benzofuran_depth = depth
                        print(f"Found benzofuran formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both transformations were found and in the correct order
    # (benzofuran formation should come after diaryl ether formation in the synthesis,
    # which means a lower depth in the retrosynthetic tree)
    if found_diaryl_ether_formation and found_benzofuran_formation:
        if benzofuran_depth < diaryl_ether_depth:
            print(
                "Strategy detected: Late-stage diaryl ether formation followed by benzofuran ring construction"
            )
            return True

    return False
