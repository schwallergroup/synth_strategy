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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis route involves an early-stage amide bond formation
    (carboxylic acid and amine to amide transformation in the early steps).
    """
    amide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found

        # Early stage is at higher depth values in retrosynthesis (depth >= 2)
        if node["type"] == "reaction" and depth >= 2:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check if this is an amide formation reaction
                    amide_formation_reaction = False

                    # Check for specific amide formation reaction types
                    if (
                        checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                        or checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        )
                        or checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction("Ester with primary amine to amide", rsmi)
                        or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                        or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                        or checker.check_reaction("Ester with ammonia to amide", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction(
                            "Acylation of secondary amines with anhydrides", rsmi
                        )
                        or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    ):
                        print(f"Amide formation reaction detected: {rsmi}")
                        amide_formation_reaction = True

                    # If not a specific reaction type, check for functional group transformation
                    if not amide_formation_reaction:
                        # Check for amide precursors in reactants
                        reactant_has_acid = any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )
                        reactant_has_acyl_halide = any(
                            checker.check_fg("Acyl halide", r) for r in reactants
                        )
                        reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                        reactant_has_anhydride = any(
                            checker.check_fg("Anhydride", r) for r in reactants
                        )

                        reactant_has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants
                        )

                        # Check for amide in product
                        product_has_amide = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        )

                        # Verify amide formation
                        if product_has_amide and (
                            (reactant_has_acid and reactant_has_amine)
                            or (reactant_has_acyl_halide and reactant_has_amine)
                            or (reactant_has_ester and reactant_has_amine)
                            or (reactant_has_anhydride and reactant_has_amine)
                        ):
                            print(f"Amide formation detected through FG analysis: {rsmi}")
                            amide_formation_reaction = True

                    if amide_formation_reaction:
                        print(f"Early-stage amide formation confirmed at depth {depth}")
                        amide_formation_found = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: early stage amide formation found = {amide_formation_found}")
    return amide_formation_found
