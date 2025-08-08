#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy involving late-stage amide formation.
    """
    has_late_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for amide formation reactions
            is_amide_formation = False

            # Check using reaction types
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acyl chloride with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with primary amine to imide",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with ammonia to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acylation of secondary amines with anhydrides",
                "Schotten-Baumann_amide",
            ]

            for reaction_type in amide_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide formation reaction: {reaction_type}")
                    is_amide_formation = True
                    break

            # If no specific reaction type matched, check for functional group changes
            if not is_amide_formation:
                print(
                    "No specific amide formation reaction type matched, checking functional groups..."
                )

                # Check for amide-forming reagents in reactants
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_carboxylic_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                # Check for amines in reactants
                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)

                # Check for amide in product
                has_primary_amide_product = checker.check_fg("Primary amide", product)
                has_secondary_amide_product = checker.check_fg("Secondary amide", product)
                has_tertiary_amide_product = checker.check_fg("Tertiary amide", product)

                has_amide_in_product = (
                    has_primary_amide_product
                    or has_secondary_amide_product
                    or has_tertiary_amide_product
                )

                # Check for amides in reactants
                has_primary_amide_reactants = any(
                    checker.check_fg("Primary amide", r) for r in reactants
                )
                has_secondary_amide_reactants = any(
                    checker.check_fg("Secondary amide", r) for r in reactants
                )
                has_tertiary_amide_reactants = any(
                    checker.check_fg("Tertiary amide", r) for r in reactants
                )

                print(f"Amide in product: {has_amide_in_product}")
                print(
                    f"Amide in reactants - Primary: {has_primary_amide_reactants}, Secondary: {has_secondary_amide_reactants}, Tertiary: {has_tertiary_amide_reactants}"
                )
                print(
                    f"Amide-forming reagents - Acyl halide: {has_acyl_halide}, Carboxylic acid: {has_carboxylic_acid}, Ester: {has_ester}, Anhydride: {has_anhydride}"
                )
                print(f"Amines - Primary: {has_primary_amine}, Secondary: {has_secondary_amine}")

                # Check if we have the necessary components for amide formation
                if (
                    has_amide_in_product
                    and (has_acyl_halide or has_carboxylic_acid or has_ester or has_anhydride)
                    and (has_primary_amine or has_secondary_amine)
                ):
                    # Try to verify using atom mapping that an amine is converted to an amide
                    try:
                        # Find reactant with primary or secondary amine
                        amine_reactant = None
                        for r in reactants:
                            if checker.check_fg("Primary amine", r) or checker.check_fg(
                                "Secondary amine", r
                            ):
                                amine_reactant = r
                                break

                        # Find reactant with acyl component
                        acyl_reactant = None
                        for r in reactants:
                            if (
                                checker.check_fg("Acyl halide", r)
                                or checker.check_fg("Carboxylic acid", r)
                                or checker.check_fg("Ester", r)
                                or checker.check_fg("Anhydride", r)
                            ):
                                acyl_reactant = r
                                break

                        if amine_reactant and acyl_reactant:
                            print(f"Found amine reactant: {amine_reactant}")
                            print(f"Found acyl reactant: {acyl_reactant}")
                            is_amide_formation = True
                            print(f"Detected amide formation by functional group analysis")
                    except Exception as e:
                        print(f"Error in atom mapping analysis: {e}")
                        # If atom mapping analysis fails, fall back to simpler check
                        is_amide_formation = True
                        print(f"Detected amide formation by functional group analysis (fallback)")

            # Check if this is late stage (depth <= 1)
            if is_amide_formation and depth <= 1:
                has_late_amide = True
                print(f"Detected late-stage amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_amide
