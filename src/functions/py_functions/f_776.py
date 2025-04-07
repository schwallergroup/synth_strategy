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
    Detects if the synthesis involves an amide formation from a carboxylic acid.
    """
    found_amide_formation = False

    def dfs_traverse(node):
        nonlocal found_amide_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide formation reaction using the checker function
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                    "Schotten-Baumann_amide",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type}")

                        # Verify that we have a carboxylic acid in reactants and amide in product
                        acid_found = False
                        amine_found = False
                        amide_found = False

                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                acid_found = True
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                amine_found = True

                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            amide_found = True

                        if acid_found and amine_found and amide_found:
                            found_amide_formation = True
                            print(
                                "Confirmed amide formation: carboxylic acid + amine → amide"
                            )
                            break

                # If no specific reaction type matched, try a more general approach
                if not found_amide_formation:
                    acid_found = False
                    amine_found = False
                    amide_found = False

                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            acid_found = True
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            amine_found = True

                    if (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        amide_found = True

                    if acid_found and amine_found and amide_found:
                        found_amide_formation = True
                        print(
                            "Found general amide formation pattern: carboxylic acid + amine → amide"
                        )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_amide_formation
