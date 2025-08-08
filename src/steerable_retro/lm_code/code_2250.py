#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects late-stage N-alkylation (in the first half of the synthesis).
    """
    n_alkylation_found = False
    max_depth = 0
    n_alkylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found, max_depth, n_alkylation_depth

        # Track maximum depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alkyl halide in reactants
            alkyl_halide_pattern = Chem.MolFromSmarts("C[Br,Cl,I,F]")

            # Check for secondary amine in reactants
            sec_amine_pattern = Chem.MolFromSmarts("[NH]")

            # Check for tertiary amine in product
            tert_amine_pattern = Chem.MolFromSmarts("[N]([C])([C])[C]")

            # Check if reactants have alkyl halide and secondary amine, and product has tertiary amine
            reactants_have_alkyl_halide = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(alkyl_halide_pattern)
                for r in reactants
                if r
            )
            reactants_have_sec_amine = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(sec_amine_pattern) for r in reactants if r
            )
            product_has_tert_amine = (
                Chem.MolFromSmiles(product).HasSubstructMatch(tert_amine_pattern)
                if product
                else False
            )

            if reactants_have_alkyl_halide and reactants_have_sec_amine and product_has_tert_amine:
                n_alkylation_found = True
                n_alkylation_depth = depth
                print(f"N-alkylation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if N-alkylation is in the first half of the synthesis (late stage)
    is_late_stage = (
        n_alkylation_found
        and n_alkylation_depth is not None
        and n_alkylation_depth <= max_depth / 2
    )

    print(
        f"Late-stage N-alkylation detected: {is_late_stage} (depth: {n_alkylation_depth}, max depth: {max_depth})"
    )
    return is_late_stage
