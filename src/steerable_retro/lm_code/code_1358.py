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
    This function detects if the synthetic route involves late-stage amide bond formation (depth < 3)
    or late-stage Boc deprotection that prepares an amine for amide coupling.
    """
    late_stage_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_found

        if node["type"] == "reaction" and depth < 3:  # Late stage = low depth
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acyl chloride with ammonia to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                # Check for Boc deprotection reactions
                boc_deprotection_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                ]

                # Check if this is an amide coupling reaction
                is_amide_coupling = False
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide coupling reaction: {reaction_type}")
                        is_amide_coupling = True
                        break

                # Check if this is a Boc deprotection reaction
                is_boc_deprotection = False
                for reaction_type in boc_deprotection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected Boc deprotection reaction: {reaction_type}")
                        is_boc_deprotection = True
                        break

                # Check for Boc group in reactants
                boc_group_in_reactants = False
                for reactant in reactants:
                    if "C(C)(C)OC(=O)N" in reactant or checker.check_fg("Carbamic ester", reactant):
                        print(f"Found Boc group in reactant: {reactant}")
                        boc_group_in_reactants = True

                # Fallback check for Boc deprotection by examining functional groups
                if not is_boc_deprotection and boc_group_in_reactants:
                    amine_in_product = checker.check_fg(
                        "Primary amine", product
                    ) or checker.check_fg("Secondary amine", product)
                    if amine_in_product:
                        print("Detected Boc deprotection by functional group analysis")
                        is_boc_deprotection = True

                # Check for amide in reactants (to verify it's newly formed)
                amide_in_reactants = False
                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                    ):
                        amide_in_reactants = True
                        print(f"Found amide in reactant: {reactant}")

                # Check for amide in product
                amide_in_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )
                if amide_in_product:
                    print(f"Found amide in product: {product}")

                # Case 1: Amide coupling reaction
                if is_amide_coupling:
                    # Verify reactants have carboxylic acid/acyl halide and amine
                    acid_found = False
                    amine_found = False

                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant) or checker.check_fg(
                            "Acyl halide", reactant
                        ):
                            print(f"Found acid/acyl halide in reactant: {reactant}")
                            acid_found = True

                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            print(f"Found amine in reactant: {reactant}")
                            amine_found = True

                    # Confirm amide is newly formed
                    if acid_found and amine_found and amide_in_product and not amide_in_reactants:
                        print(f"Late-stage amide coupling confirmed at depth {depth}")
                        late_stage_amide_found = True

                # Case 2: Boc deprotection (preparing amine for amide coupling)
                elif is_boc_deprotection:
                    # Check if product has a primary or secondary amine
                    if checker.check_fg("Primary amine", product) or checker.check_fg(
                        "Secondary amine", product
                    ):
                        print(f"Late-stage Boc deprotection revealing amine at depth {depth}")

                        # If the product already contains amide groups, this is part of an amide strategy
                        if amide_in_product:
                            print(
                                f"Product contains amide groups, confirming relevance to amide strategy"
                            )
                            late_stage_amide_found = True
                        else:
                            # Check if this amine will be used in a subsequent amide coupling
                            for child in node.get("children", []):
                                if child["type"] == "mol" and (
                                    checker.check_fg("Primary amine", child["smiles"])
                                    or checker.check_fg("Secondary amine", child["smiles"])
                                ):
                                    # This is a deprotected amine that could be used in amide coupling
                                    print(
                                        f"Deprotected amine could be used in subsequent amide coupling"
                                    )
                                    late_stage_amide_found = True
                                    break

                # Case 3: Alternative check if reaction type check failed but we have acid/acyl halide and amine
                elif not is_amide_coupling and not is_boc_deprotection:
                    acid_found = False
                    amine_found = False

                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant) or checker.check_fg(
                            "Acyl halide", reactant
                        ):
                            print(f"Found acid/acyl halide in reactant: {reactant}")
                            acid_found = True

                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            print(f"Found amine in reactant: {reactant}")
                            amine_found = True

                    # Confirm amide is newly formed
                    if acid_found and amine_found and amide_in_product and not amide_in_reactants:
                        print(f"Late-stage amide formation detected at depth {depth}")
                        late_stage_amide_found = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Late-stage amide coupling strategy found: {late_stage_amide_found}")

    return late_stage_amide_found
