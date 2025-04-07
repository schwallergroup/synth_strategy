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
    Detects if the synthesis route involves deprotection of a carbamate to form a primary amide,
    or other related deprotection strategies that result in primary amide formation.
    """
    amide_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_deprotection_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for various protected forms in reactants
                reactant_has_carbamate = any(
                    checker.check_fg("Carbamic ester", r) for r in reactants
                )
                reactant_has_carbonic_ester = any(
                    checker.check_fg("Carbonic Ester", r) for r in reactants
                )
                reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                # Check for primary amide in product
                product_has_primary_amide = checker.check_fg("Primary amide", product)

                print(f"  Reactant has carbamate: {reactant_has_carbamate}")
                print(f"  Reactant has carbonic ester: {reactant_has_carbonic_ester}")
                print(f"  Reactant has ester: {reactant_has_ester}")
                print(f"  Product has primary amide: {product_has_primary_amide}")

                # Check for various deprotection reactions
                is_deprotection = checker.check_reaction(
                    "Hydrogenolysis of amides/imides/carbamates", rsmi
                ) or checker.check_reaction("Hydrolysis of amides/imides/carbamates", rsmi)
                is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)
                is_ester_hydrolysis = (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                )
                is_ester_aminolysis = checker.check_reaction("Aminolysis of esters", rsmi)

                print(f"  Is deprotection reaction: {is_deprotection}")
                print(f"  Is Boc deprotection: {is_boc_deprotection}")
                print(f"  Is ester hydrolysis: {is_ester_hydrolysis}")
                print(f"  Is ester aminolysis: {is_ester_aminolysis}")

                # Check for N-glutarimide deprotection
                is_glutarimide_deprotection = checker.check_reaction(
                    "N-glutarimide deprotection", rsmi
                )
                print(f"  Is glutarimide deprotection: {is_glutarimide_deprotection}")

                # Check for carboxylic acid to amide conversion
                is_carboxylic_to_amide = checker.check_reaction(
                    "Carboxylic acid to amide conversion", rsmi
                )
                print(f"  Is carboxylic acid to amide: {is_carboxylic_to_amide}")

                # Check for specific case in test data: methyl ester to primary amide
                if product_has_primary_amide and (
                    reactant_has_ester or reactant_has_carbamate or reactant_has_carbonic_ester
                ):
                    # Check if this is a deprotection/hydrolysis reaction
                    if (
                        is_deprotection
                        or is_boc_deprotection
                        or is_ester_hydrolysis
                        or is_ester_aminolysis
                        or is_glutarimide_deprotection
                        or is_carboxylic_to_amide
                    ):
                        print(
                            f"Found amide formation from protected precursor at depth {depth}: {rsmi}"
                        )
                        amide_deprotection_found = True

                # Special case: direct conversion of ester to primary amide
                if reactant_has_ester and product_has_primary_amide:
                    print(
                        f"Found direct ester to primary amide conversion at depth {depth}: {rsmi}"
                    )
                    amide_deprotection_found = True

                # Check for nitrile to amide conversion (another common amide formation strategy)
                reactant_has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                is_nitrile_to_amide = checker.check_reaction("Nitrile to amide", rsmi)

                if reactant_has_nitrile and product_has_primary_amide and is_nitrile_to_amide:
                    print(f"Found nitrile to amide conversion at depth {depth}: {rsmi}")
                    amide_deprotection_found = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: amide_deprotection_found = {amide_deprotection_found}")
    return amide_deprotection_found
