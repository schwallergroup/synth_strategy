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
    This function detects if the synthesis includes a late-stage functional group interconversion,
    specifically focusing on transformations in the final steps.
    """
    late_stage_detected = False

    # List of functional group interconversions to check
    fg_interconversions = [
        # Format: (starting_fg, ending_fg, reaction_type)
        ("Aldehyde", "Alkyne", ""),  # Corey-Fuchs reaction not in provided list
        ("Carboxylic acid", "Ester", "Esterification of Carboxylic Acids"),
        ("Primary alcohol", "Primary halide", "Alkyl chlorides from alcohols"),
        ("Nitrile", "Primary amine", "Reduction of nitrile to amine"),
        ("Ketone", "Secondary alcohol", "Reduction of ketone to secondary alcohol"),
        ("Ester", "Primary alcohol", "Reduction of ester to primary alcohol"),
        ("Nitro group", "Primary amine", "Reduction of nitro groups to amines"),
        (
            "Primary alcohol",
            "Aldehyde",
            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        ),
        (
            "Secondary alcohol",
            "Ketone",
            "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        ),
        ("Alkyne", "Ketone", "Hydration of alkyne to ketone"),
        ("Alkyne", "Aldehyde", "Hydration of alkyne to aldehyde"),
        ("Primary halide", "Primary alcohol", "Primary alkyl halide to alcohol"),
        ("Secondary halide", "Secondary alcohol", "Secondary alkyl halide to alcohol"),
        ("Aldehyde", "Primary alcohol", "Reduction of aldehydes and ketones to alcohols"),
        ("Ketone", "Secondary alcohol", "Reduction of aldehydes and ketones to alcohols"),
        ("Carboxylic acid", "Primary alcohol", "Reduction of carboxylic acid to primary alcohol"),
        (
            "Ester",
            "Carboxylic acid",
            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        ),
        ("Nitrile", "Carboxylic acid", "Oxidation of nitrile to carboxylic acid"),
        ("Nitrile", "Amide", "Nitrile to amide"),
        ("Alcohol", "Ether", "Alcohol to ether"),
        ("Alcohol", "Azide", "Alcohol to azide"),
        ("Primary amine", "Azide", "Amine to azide"),
        ("Carboxylic acid", "Amide", "Carboxylic acid with primary amine to amide"),
        ("Ester", "Amide", "Ester with primary amine to amide"),
        ("Acyl halide", "Amide", "Acyl chloride with primary amine to amide (Schotten-Baumann)"),
        ("Acyl halide", "Ester", "Schotten-Baumann to ester"),
        ("Alcohol", "Triflate", "Alcohol to triflate conversion"),
        ("Alcohol", "Tosylate", "Formation of Sulfonic Esters"),
        ("Alcohol", "Mesylate", "Formation of Sulfonic Esters"),
        ("Boronic acid", "Phenol", "Oxidation of boronic acids"),
        ("Boronic ester", "Phenol", "Oxidation of boronic esters"),
        ("Aromatic halide", "Boronic acid", "Preparation of boronic acids"),
        ("Aromatic halide", "Boronic ester", "Preparation of boronic ethers"),
    ]

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_detected

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check for functional group interconversions
            for start_fg, end_fg, reaction_name in fg_interconversions:
                # Find a reactant with the starting functional group
                reactant_with_fg = None
                for reactant in reactants:
                    if checker.check_fg(start_fg, reactant):
                        reactant_with_fg = reactant
                        break

                # Check if product has the ending functional group
                has_end_fg = checker.check_fg(end_fg, product)

                # Only proceed if we found a reactant with the starting functional group
                has_start_fg = reactant_with_fg is not None

                # Check if this is the expected reaction type (if provided)
                is_expected_reaction = True
                if reaction_name:
                    is_expected_reaction = checker.check_reaction(reaction_name, rsmi)

                # If we have the starting FG, ending FG, and correct reaction type (if specified)
                if has_start_fg and has_end_fg and is_expected_reaction:
                    print(f"Found potential {start_fg} to {end_fg} transformation")
                    print(f"Reaction type check: {is_expected_reaction}")

                    # If depth is 0, 1, or 2, it's considered late-stage
                    if depth <= 2:
                        print(
                            f"Detected late-stage {start_fg} to {end_fg} transformation at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        late_stage_detected = True
                    else:
                        print(
                            f"Found {start_fg} to {end_fg} transformation but at depth {depth} (not late-stage)"
                        )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_detected
