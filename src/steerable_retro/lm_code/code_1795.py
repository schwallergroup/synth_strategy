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
    This function detects a synthetic strategy involving multiple amide bond disconnections
    in the retrosynthetic route.
    """
    # Count amide disconnections
    amide_disconnection_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_disconnection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reaction_id = node.get("metadata", {}).get("ID", "unknown")

            # Check for amide formation/disconnection reactions using the checker function
            is_amide_reaction = False
            matched_reaction_type = None

            # Check for specific amide formation/disconnection reaction types
            amide_reaction_types = [
                # Amide formation reactions (appear as disconnections in retrosynthesis)
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Ester with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Acyl chloride with ammonia to amide",
                "Schotten-Baumann_amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acylation of secondary amines with anhydrides",
                "Carboxylic acid to amide conversion",
                # Amide hydrolysis/cleavage reactions (direct disconnections)
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Hydrogenolysis of amides/imides/carbamates",
                "Hydrolysis of amides/imides/carbamates",
            ]

            for reaction_type in amide_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    is_amide_reaction = True
                    matched_reaction_type = reaction_type
                    break

            # If no specific reaction type matched, check for functional group patterns
            if not is_amide_reaction:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                product = product_part

                # Check for amide formation (in retrosynthesis: product has acid/amine, reactants have amide)
                has_amide_in_reactants = any(
                    checker.check_fg("Primary amide", reactant)
                    or checker.check_fg("Secondary amide", reactant)
                    or checker.check_fg("Tertiary amide", reactant)
                    for reactant in reactants
                )

                has_acid_in_product = checker.check_fg("Carboxylic acid", product)
                has_amine_in_product = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                # Check for amide formation (in forward direction: reactants have acid/amine, product has amide)
                has_amide_in_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                has_acid_in_reactants = any(
                    checker.check_fg("Carboxylic acid", reactant) for reactant in reactants
                )
                has_amine_in_reactants = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    or checker.check_fg("Tertiary amine", reactant)
                    for reactant in reactants
                )

                # Either direction could indicate an amide disconnection in retrosynthesis
                if (has_amide_in_reactants and (has_acid_in_product or has_amine_in_product)) or (
                    has_amide_in_product and has_acid_in_reactants and has_amine_in_reactants
                ):
                    is_amide_reaction = True
                    matched_reaction_type = "Custom amide pattern"

            if is_amide_reaction:
                amide_disconnection_count += 1
                print(
                    f"Found amide disconnection at ID {reaction_id} - Reaction type: {matched_reaction_type}"
                )
                print(f"Reaction SMILES: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Total amide disconnections found: {amide_disconnection_count}")

    # Strategy requires at least 2 amide disconnections
    return amide_disconnection_count >= 2
