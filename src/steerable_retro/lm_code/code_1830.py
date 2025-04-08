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
    This function detects orthogonal protection of carboxylic acids using t-butyl and methyl esters
    with preservation of an azide functional group.
    """
    # Track if we've found the key features
    has_tbutyl_protection = False
    has_methyl_protection = False
    has_azide_preservation = False
    has_amide_bond = False
    has_deprotection = False

    # Track molecules at different depths to check for preservation
    molecules_with_azide = set()

    # Track molecules with protected carboxylic acids
    molecules_with_tbutyl_ester = set()
    molecules_with_methyl_ester = set()

    def dfs_traverse(node, depth=0):
        nonlocal has_tbutyl_protection, has_methyl_protection, has_azide_preservation
        nonlocal has_amide_bond, has_deprotection, molecules_with_azide
        nonlocal molecules_with_tbutyl_ester, molecules_with_methyl_ester

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for azide group in molecules
            if checker.check_fg("Azide", mol_smiles):
                molecules_with_azide.add(depth)
                print(f"Found azide at depth {depth}: {mol_smiles}")

            # Check for t-butyl ester in molecules
            if checker.check_fg("Ester", mol_smiles) and "OC(C)(C)C" in mol_smiles:
                molecules_with_tbutyl_ester.add(depth)
                print(f"Found t-butyl ester at depth {depth}: {mol_smiles}")

            # Check for methyl ester in molecules
            if checker.check_fg("Ester", mol_smiles) and "COC(=O)" in mol_smiles:
                molecules_with_methyl_ester.add(depth)
                print(f"Found methyl ester at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check for t-butyl ester protection
                if (
                    checker.check_fg("Carboxylic acid", reactants_part)
                    and "OC(C)(C)C" in product
                    and not any("OC(C)(C)C" in r for r in reactants)
                ):
                    has_tbutyl_protection = True
                    print(f"Found t-butyl ester protection at depth {depth}: {rsmi}")

                # Check for methyl ester protection
                if (
                    checker.check_fg("Carboxylic acid", reactants_part)
                    and "COC(=O)" in product
                    and not any("COC(=O)" in r for r in reactants)
                ):
                    has_methyl_protection = True
                    print(f"Found methyl ester protection at depth {depth}: {rsmi}")

                # Check for protection of carboxylic acid (general)
                if checker.check_reaction("Protection of carboxylic acid", rsmi):
                    if "OC(C)(C)C" in product and not any("OC(C)(C)C" in r for r in reactants):
                        has_tbutyl_protection = True
                        print(
                            f"Found t-butyl protection via reaction check at depth {depth}: {rsmi}"
                        )
                    if "COC(=O)" in product and not any("COC(=O)" in r for r in reactants):
                        has_methyl_protection = True
                        print(
                            f"Found methyl protection via reaction check at depth {depth}: {rsmi}"
                        )

                # Check for amide bond formation
                if (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                ):
                    if not any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants
                    ):
                        has_amide_bond = True
                        print(f"Found amide bond formation at depth {depth}: {rsmi}")

                # Check for ester deprotection
                if checker.check_fg("Carboxylic acid", product):
                    if any(checker.check_fg("Ester", r) for r in reactants):
                        if (
                            checker.check_reaction("COOH ethyl deprotection", rsmi)
                            or checker.check_reaction(
                                "Ester saponification (alkyl deprotection)", rsmi
                            )
                            or checker.check_reaction(
                                "Ester saponification (methyl deprotection)", rsmi
                            )
                            or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                            or checker.check_reaction(
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                rsmi,
                            )
                        ):
                            has_deprotection = True
                            print(f"Found ester deprotection at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if azide is preserved across multiple steps
    if len(molecules_with_azide) >= 2:
        has_azide_preservation = True
        print(f"Azide group preserved across depths: {molecules_with_azide}")

    # Check if we have both types of protection in different molecules
    orthogonal_protection = (
        len(molecules_with_tbutyl_ester) > 0 and len(molecules_with_methyl_ester) > 0
    )
    if orthogonal_protection:
        print(
            f"Found orthogonal protection: t-butyl at depths {molecules_with_tbutyl_ester}, methyl at depths {molecules_with_methyl_ester}"
        )

    # Return True if all key features are present
    # We need orthogonal protection (both t-butyl and methyl esters) and azide preservation
    strategy_present = orthogonal_protection and has_azide_preservation

    print(f"Orthogonal carboxylic acid protection strategy detected: {strategy_present}")
    print(f"t-butyl protection: {has_tbutyl_protection}")
    print(f"methyl protection: {has_methyl_protection}")
    print(f"azide preservation: {has_azide_preservation}")
    print(f"deprotection: {has_deprotection}")
    print(f"amide bond: {has_amide_bond}")
    print(f"orthogonal protection: {orthogonal_protection}")

    return strategy_present
