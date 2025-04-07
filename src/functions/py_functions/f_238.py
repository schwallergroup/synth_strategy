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
    Detects if the synthesis route contains an amide formation from a primary amine.
    This includes carbamate formation (like Boc protection).
    """
    amide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation reactions directly
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Acylation of primary amines",
                "Boc amine protection",
                "Boc amine protection explicit",
                "Boc amine protection with Boc anhydride",
                "Boc amine protection (ethyl Boc)",
                "Boc amine protection of primary amine",
            ]

            for reaction_name in amide_formation_reactions:
                if checker.check_reaction(reaction_name, rsmi):
                    # Verify that a primary amine is involved
                    if any(checker.check_fg("Primary amine", r) for r in reactants):
                        # Verify that an amide or carbamate is formed
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                            or checker.check_fg("Carbamic ester", product)
                            or checker.check_fg("Carbamic acid", product)
                        ):
                            amide_formation_found = True
                            print(f"Amide formation found at depth {depth}")
                            print(f"Reaction: {reaction_name}")
                            print(f"RSMI: {rsmi}")
                            break

            # If no specific reaction was found, check for general pattern
            if not amide_formation_found:
                # Check if reactants contain primary amine
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )

                # Check if product contains amide or carbamate
                has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                    or checker.check_fg("Carbamic ester", product)
                    or checker.check_fg("Carbamic acid", product)
                )

                if has_primary_amine and has_amide:
                    print(f"Potential amide formation detected at depth {depth}")
                    print(f"RSMI: {rsmi}")

                    # Check for acyl halides, carboxylic acids, esters, or Boc sources in reactants
                    has_acyl_source = any(
                        checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Ester", r)
                        or checker.check_fg("Anhydride", r)
                        or "BOC" in r.upper()
                        or "OC(=O)O" in r
                        for r in reactants
                    )

                    if has_acyl_source:
                        amide_formation_found = True
                        print("Confirmed as amide formation from primary amine")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return amide_formation_found
