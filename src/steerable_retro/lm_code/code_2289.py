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
    Detects if the synthesis uses multiple functional group interconversions
    (e.g., ester to alcohol, amide to amine).
    """
    interconversions = 0
    interconversion_types = set()  # Track unique types of interconversions

    def dfs_traverse(node):
        nonlocal interconversions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for various functional group interconversions
            found_interconversion = False

            # 1. Ester to alcohol reduction
            if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                print(f"Found ester to primary alcohol reduction: {rsmi}")
                interconversion_types.add("ester_to_alcohol")
                found_interconversion = True

            # 2. Amide to amine reduction
            elif any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                ]
            ):
                print(f"Found amide to amine reduction: {rsmi}")
                interconversion_types.add("amide_to_amine")
                found_interconversion = True

            # 3. Nitrile to amine reduction
            elif checker.check_reaction("Reduction of nitrile to amine", rsmi):
                print(f"Found nitrile to amine reduction: {rsmi}")
                interconversion_types.add("nitrile_to_amine")
                found_interconversion = True

            # 4. Carboxylic acid to ester
            elif checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                print(f"Found carboxylic acid to ester: {rsmi}")
                interconversion_types.add("acid_to_ester")
                found_interconversion = True

            # 5. Alcohol to ether
            elif checker.check_reaction("Alcohol to ether", rsmi) or checker.check_reaction(
                "Williamson Ether Synthesis", rsmi
            ):
                print(f"Found alcohol to ether: {rsmi}")
                interconversion_types.add("alcohol_to_ether")
                found_interconversion = True

            # 6. Alcohol to halide
            elif any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Alcohol to chloride_SOCl2",
                    "Primary alkyl halide to alcohol",
                    "Secondary alkyl halide to alcohol",
                    "Alkyl chlorides from alcohols",
                    "Alkyl bromides from alcohols",
                    "Alkyl iodides from alcohols",
                ]
            ):
                print(f"Found alcohol to halide conversion: {rsmi}")
                interconversion_types.add("alcohol_to_halide")
                found_interconversion = True

            # 7. Aldehyde/Ketone to alcohol
            elif checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi):
                print(f"Found carbonyl to alcohol reduction: {rsmi}")
                interconversion_types.add("carbonyl_to_alcohol")
                found_interconversion = True

            # 8. Nitro to amine reduction
            elif checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Found nitro to amine reduction: {rsmi}")
                interconversion_types.add("nitro_to_amine")
                found_interconversion = True

            # 9. Alcohol to carboxylic acid oxidation
            elif checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi):
                print(f"Found alcohol to carboxylic acid oxidation: {rsmi}")
                interconversion_types.add("alcohol_to_acid")
                found_interconversion = True

            # 10. Amine to amide
            elif any(
                checker.check_reaction(rxn, rsmi)
                for rxn in [
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Carboxylic acid with primary amine to amide",
                ]
            ):
                print(f"Found amine to amide conversion: {rsmi}")
                interconversion_types.add("amine_to_amide")
                found_interconversion = True

            # If no specific reaction check matched, try checking for functional group changes
            if not found_interconversion:
                for reactant in reactants:
                    # Check for functional group changes by comparing reactant and product
                    if checker.check_fg("Ester", reactant) and (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    ):
                        print(f"Found ester to alcohol conversion (FG check): {rsmi}")
                        interconversion_types.add("ester_to_alcohol")
                        found_interconversion = True
                        break

                    elif (
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                    ) and (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    ):
                        print(f"Found amide to amine conversion (FG check): {rsmi}")
                        interconversion_types.add("amide_to_amine")
                        found_interconversion = True
                        break

                    elif checker.check_fg("Nitrile", reactant) and (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    ):
                        print(f"Found nitrile to amine conversion (FG check): {rsmi}")
                        interconversion_types.add("nitrile_to_amine")
                        found_interconversion = True
                        break

                    elif checker.check_fg("Carboxylic acid", reactant) and checker.check_fg(
                        "Ester", product
                    ):
                        print(f"Found carboxylic acid to ester conversion (FG check): {rsmi}")
                        interconversion_types.add("acid_to_ester")
                        found_interconversion = True
                        break

            if found_interconversion:
                interconversions += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(
        f"Found {interconversions} functional group interconversions of {len(interconversion_types)} different types"
    )
    # Return True if at least 2 interconversions are found
    return interconversions >= 2
