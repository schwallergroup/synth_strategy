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
    Detects a synthesis featuring multiple nitrogen protecting group manipulations
    (≥2 different protecting groups) with a late-stage introduction of a heteroaromatic ring
    via nucleophilic aromatic substitution.
    """
    # Track protecting groups and late-stage heteroaryl introduction
    protecting_groups_used = set()
    late_stage_heteroaryl_introduction = False
    max_depth = 0

    # SMARTS patterns for protecting groups
    phthalimide_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1C(=O)[N]C(=O)")
    boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])OC(=O)[N]")
    cbz_pattern = Chem.MolFromSmarts("O=C([O][C][c]1[cH][cH][cH][cH][cH]1)[N]")

    # Pattern for heteroaromatic rings
    heteroaromatic_pattern = Chem.MolFromSmarts("[a;!c]1[a][a][a][a][a]1")

    def dfs_traverse(node, depth=0):
        nonlocal protecting_groups_used, late_stage_heteroaryl_introduction, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for protecting groups
                if mol.HasSubstructMatch(phthalimide_pattern):
                    protecting_groups_used.add("phthalimide")
                if mol.HasSubstructMatch(boc_pattern):
                    protecting_groups_used.add("boc")
                if mol.HasSubstructMatch(cbz_pattern):
                    protecting_groups_used.add("cbz")

        elif node["type"] == "reaction":
            # Check for late-stage heteroaryl introduction
            if depth <= 1:  # Late stage (depth 0 or 1)
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has a heteroaromatic ring
                    reactant_has_heteroaromatic = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(heteroaromatic_pattern):
                            reactant_has_heteroaromatic = True

                    # Check if product has a heteroaromatic ring
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(
                        heteroaromatic_pattern
                    ):
                        if reactant_has_heteroaromatic:
                            # Check if this is likely an SNAr reaction (halogen on heteroaromatic)
                            halogen_heteroaromatic = Chem.MolFromSmarts(
                                "[a;!c][a]([F,Cl,Br,I])"
                            )
                            for reactant in reactants:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol and mol.HasSubstructMatch(
                                    halogen_heteroaromatic
                                ):
                                    late_stage_heteroaryl_introduction = True
                                    print(
                                        f"Late-stage heteroaryl introduction detected at depth {depth}"
                                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    has_multiple_protecting_groups = len(protecting_groups_used) >= 2

    if has_multiple_protecting_groups and late_stage_heteroaryl_introduction:
        print(
            f"Protection-deprotection cascade with late-stage heteroaryl introduction detected"
        )
        print(f"Protecting groups used: {protecting_groups_used}")
        return True
    return False
