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
    This function detects multiple amide coupling reactions in the synthesis.
    """
    amide_couplings = 0

    def dfs_traverse(node):
        nonlocal amide_couplings

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check for amide coupling reactions using the checker function
            amide_coupling_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acyl chloride with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with ammonia to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Carboxylic acid to amide conversion",
            ]

            for reaction_type in amide_coupling_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found amide coupling reaction: {reaction_type}")
                    print(f"Reaction SMILES: {rsmi}")
                    amide_couplings += 1
                    break

            # If no specific reaction type matched, check for functional groups
            if not any(
                checker.check_reaction(reaction_type, rsmi)
                for reaction_type in amide_coupling_reactions
            ):
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                has_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                )

                # Check for amines in reactants
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )
                has_amine = has_primary_amine or has_secondary_amine

                # Check for amide in product
                has_primary_amide = checker.check_fg("Primary amide", product)
                has_secondary_amide = checker.check_fg("Secondary amide", product)
                has_tertiary_amide = checker.check_fg("Tertiary amide", product)
                has_amide = (
                    has_primary_amide or has_secondary_amide or has_tertiary_amide
                )

                if has_acid and has_amine and has_amide:
                    print(f"Found amide coupling reaction by functional group analysis")
                    print(f"Reaction SMILES: {rsmi}")
                    amide_couplings += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total amide coupling reactions found: {amide_couplings}")

    return amide_couplings >= 2
