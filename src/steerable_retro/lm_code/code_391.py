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
    This function detects if the synthesis route involves amide bond formation.
    """
    amide_formation = False

    def dfs_traverse(node):
        nonlocal amide_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for specific amide formation reaction types
                amide_formation_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Acyl chloride with ammonia to amide",
                    "Schotten-Baumann_amide",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        amide_formation = True
                        print(f"Amide formation detected via {reaction_type}: {rsmi}")
                        return

                # Fallback to pattern matching if specific reaction types aren't detected
                # Check for carboxylic acid in reactants
                has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)

                # Check for acyl halide in reactants (another common amide precursor)
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)

                # Check for ester in reactants (can form amides with amines)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                # Check for amine in reactants
                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)
                has_amine = has_primary_amine or has_secondary_amine

                # Check for amide in product
                has_primary_amide = checker.check_fg("Primary amide", product)
                has_secondary_amide = checker.check_fg("Secondary amide", product)
                has_tertiary_amide = checker.check_fg("Tertiary amide", product)
                has_amide_product = has_primary_amide or has_secondary_amide or has_tertiary_amide

                # Check for amide formation conditions
                if (has_acid or has_acyl_halide or has_ester) and has_amine and has_amide_product:
                    amide_formation = True
                    print(f"Amide formation detected through pattern matching: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    return amide_formation
