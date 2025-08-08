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
    This function detects amide bond formation from carboxylic acid and amine.
    """
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check for required functional groups in reactants
                has_carboxylic_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants if r
                )
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants if r
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants if r
                )
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants if r)
                has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants if r)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants if r)

                print(
                    f"Reactant functional groups: Carboxylic acid: {has_carboxylic_acid}, Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}, Acyl halide: {has_acyl_halide}, Anhydride: {has_anhydride}, Ester: {has_ester}"
                )

                # Check for amide in product
                has_primary_amide = checker.check_fg("Primary amide", product)
                has_secondary_amide = checker.check_fg("Secondary amide", product)
                has_tertiary_amide = checker.check_fg("Tertiary amide", product)
                has_amide = has_primary_amide or has_secondary_amide or has_tertiary_amide

                print(
                    f"Product functional groups: Primary amide: {has_primary_amide}, Secondary amide: {has_secondary_amide}, Tertiary amide: {has_tertiary_amide}"
                )

                # Check for specific amide formation reactions
                is_amide_formation_reaction = (
                    checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                )

                print(f"Is amide formation reaction: {is_amide_formation_reaction}")

                # If we have the right functional groups and reaction type
                if (
                    has_amide
                    and (
                        (has_carboxylic_acid and (has_primary_amine or has_secondary_amine))
                        or (has_acyl_halide and (has_primary_amine or has_secondary_amine))
                        or (has_anhydride and (has_primary_amine or has_secondary_amine))
                        or (has_ester and (has_primary_amine or has_secondary_amine))
                    )
                    and is_amide_formation_reaction
                ):
                    has_amide_formation = True
                    print(f"Detected amide bond formation: {rsmi}")

                # Check for amide formation even if specific reaction type not detected
                elif has_amide and (
                    (has_carboxylic_acid and (has_primary_amine or has_secondary_amine))
                    or (has_acyl_halide and (has_primary_amine or has_secondary_amine))
                    or (has_anhydride and (has_primary_amine or has_secondary_amine))
                    or (has_ester and (has_primary_amine or has_secondary_amine))
                ):
                    # If we have the right functional groups but reaction type not detected,
                    # check if product has amide and reactants don't
                    reactants_have_amide = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants
                        if r
                    )

                    if has_amide and not reactants_have_amide:
                        has_amide_formation = True
                        print(
                            f"Detected amide bond formation (by functional group analysis): {rsmi}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return has_amide_formation
