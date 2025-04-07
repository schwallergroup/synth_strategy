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
    This function detects if N-acylation is used as a late-stage functionalization.
    It looks for an amide formation in the last two steps of synthesis.
    """
    n_acylation_at_depth_0_or_1 = False

    def dfs_traverse(node, depth=0):
        nonlocal n_acylation_at_depth_0_or_1

        if node["type"] == "reaction" and depth <= 1:  # Only check the last two steps
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for N-acylation reactions directly
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                ):

                    # Verify that we have the right functional groups in reactants and products
                    has_amine = False
                    has_acylating_agent = False
                    has_amide_product = False

                    # Check reactants for amines and acylating agents
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine = True

                        if (
                            checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_acylating_agent = True

                    # Check product for amide
                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        has_amide_product = True

                    if has_amine and has_acylating_agent and has_amide_product:
                        print(f"Found N-acylation at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        n_acylation_at_depth_0_or_1 = True

                # Fallback check for N-acylation if reaction type check fails
                elif not n_acylation_at_depth_0_or_1:
                    has_amine = False
                    has_acylating_agent = False
                    has_amide_product = False

                    # Check reactants for amines and acylating agents
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine = True

                        if (
                            checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_acylating_agent = True

                    # Check product for amide
                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        has_amide_product = True

                    if has_amine and has_acylating_agent and has_amide_product:
                        print(f"Found N-acylation (fallback check) at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        n_acylation_at_depth_0_or_1 = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return n_acylation_at_depth_0_or_1
