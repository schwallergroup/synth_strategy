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
    This function detects if the synthesis involves late-stage incorporation of a pyrrolidine ring.
    """
    late_stage_pyrrolidine = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_pyrrolidine

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Consider only reactions in the first 3 levels (late stage)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains pyrrolidine
                has_pyrrolidine_product = checker.check_ring("pyrrolidine", product)

                # Check if any reactant contains pyrrolidine
                pyrrolidine_reactants = [
                    r for r in reactants if checker.check_ring("pyrrolidine", r)
                ]

                # Check for pyrrolidine-forming reactions
                pyrrolidine_forming_reaction = (
                    checker.check_reaction("Paal-Knorr pyrrole synthesis", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("reductive amination with ketone", rsmi)
                )

                # Check for direct incorporation of pyrrolidine
                pyrrolidine_incorporation = (
                    has_pyrrolidine_product and len(pyrrolidine_reactants) > 0
                )

                # Check if the pyrrolidine is actually being incorporated (not just present in both)
                if pyrrolidine_incorporation or (
                    has_pyrrolidine_product and pyrrolidine_forming_reaction
                ):
                    print(
                        f"Detected late-stage pyrrolidine incorporation at depth {depth}"
                    )
                    print(f"Reaction SMILES: {rsmi}")
                    late_stage_pyrrolidine = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_pyrrolidine
