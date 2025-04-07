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
    Detects if the synthesis follows a convergent approach with a late-stage amide coupling.
    """
    convergent_synthesis = False
    late_stage_amide = False

    # Define amide coupling reaction types to check
    amide_coupling_reactions = [
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Carboxylic acid with primary amine to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Schotten-Baumann_amide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis, late_stage_amide

        if node["type"] == "reaction" and depth <= 1:  # Check top two levels (late-stage)
            # Check if reaction has multiple reactants (convergent)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                if len(reactants) >= 2:
                    print(f"Found convergent step at depth {depth} with {len(reactants)} reactants")
                    convergent_synthesis = True

                    # Check for amide coupling reaction
                    for reaction_type in amide_coupling_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Found amide coupling reaction: {reaction_type}")
                            late_stage_amide = True
                            break

                    # If reaction type check failed, check for functional groups
                    if not late_stage_amide:
                        # Check for carboxylic acid in reactants
                        has_acid = any(
                            checker.check_fg("Carboxylic acid", reactant) for reactant in reactants
                        )

                        # Check for amine in reactants
                        has_amine = any(
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            for reactant in reactants
                        )

                        # Check for amide in product
                        has_amide = (
                            checker.check_fg("Primary amide", products_part)
                            or checker.check_fg("Secondary amide", products_part)
                            or checker.check_fg("Tertiary amide", products_part)
                        )

                        if has_acid and has_amine and has_amide:
                            print(
                                f"Found amide coupling based on functional groups at depth {depth}"
                            )
                            late_stage_amide = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Convergent synthesis: {convergent_synthesis}")
    print(f"Late-stage amide coupling: {late_stage_amide}")

    return convergent_synthesis and late_stage_amide
