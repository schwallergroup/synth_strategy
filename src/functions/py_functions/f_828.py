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
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects a synthetic strategy involving SNAr reactions with piperazine
    as a nucleophile.
    """
    has_snar_with_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar_with_piperazine

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an SNAr reaction
                is_snar = (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                )

                # Check for aromatic halide in reactants
                has_aromatic_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                )

                # Check for piperazine in reactants
                has_piperazine = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )

                # Check for aryl-piperazine in product
                has_aryl_piperazine = checker.check_ring("piperazine", product)

                if has_aromatic_halide and has_piperazine and has_aryl_piperazine:
                    print(f"Detected SNAr with piperazine at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    has_snar_with_piperazine = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_snar_with_piperazine
