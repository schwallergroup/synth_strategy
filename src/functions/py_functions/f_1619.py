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
    This function detects if the route follows a "build-couple-pair" strategy,
    where a complex scaffold is built through sequential coupling reactions,
    with functional group manipulations in late stages.
    """
    # We need to check for:
    # 1. At least two coupling reactions (e.g., Suzuki)
    # 2. Late-stage functional group manipulations
    # 3. Increasing ring count throughout the synthesis

    coupling_reactions = 0
    late_stage_fg_manipulation = False

    def dfs_traverse(node):
        nonlocal coupling_reactions, late_stage_fg_manipulation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            depth = node.get("metadata", {}).get("depth", -1)

            # Check for coupling reaction (simplified)
            has_boron = False
            has_halide = False

            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                        ) or mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6]-[#5](-[#8;R])-[#8;R]")
                        ):
                            has_boron = True

                        if mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[#6;a]-[#35,#53,#17]")
                        ):
                            has_halide = True

            if has_boron and has_halide:
                coupling_reactions += 1
                print(f"Found coupling reaction at depth {depth}: {rsmi}")

            # Check for late-stage functional group manipulation
            if depth <= 1:  # Late stage
                # Check for functional group changes (simplified)
                if any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]-[#9]")
                    )
                    for r in reactants
                    if r
                ) and Chem.MolFromSmiles(product).HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                ):
                    late_stage_fg_manipulation = True
                    print(f"Found late-stage FG manipulation: {rsmi}")

                # Check for deprotection
                if any(
                    Chem.MolFromSmiles(r).HasSubstructMatch(
                        Chem.MolFromSmarts("[#14]([#6])([#6])([#6])[#8]-[#6]")
                    )
                    for r in reactants
                    if r
                ) and Chem.MolFromSmiles(product).HasSubstructMatch(
                    Chem.MolFromSmarts("[#8;H1]-[#6]")
                ):
                    late_stage_fg_manipulation = True
                    print(f"Found late-stage deprotection: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    build_couple_pair = coupling_reactions >= 2 and late_stage_fg_manipulation
    print(
        f"Build-couple-pair strategy: {build_couple_pair} (Couplings: {coupling_reactions}, Late-stage FG: {late_stage_fg_manipulation})"
    )
    return build_couple_pair
