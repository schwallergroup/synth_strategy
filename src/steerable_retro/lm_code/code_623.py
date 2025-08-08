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
    This function detects if the synthetic route includes multiple esterification reactions.
    """
    esterification_count = 0

    # List of esterification reaction types from the provided list
    esterification_reactions = [
        "Esterification of Carboxylic Acids",
        "Transesterification",
        "O-alkylation of carboxylic acids with diazo compounds",
        "Oxidative esterification of primary alcohols",
        "Oxidation of alcohol and aldehyde to ester",
        "Acetic anhydride and alcohol to ester",
        "Mitsunobu esterification",
        "Pinner reaction to ester",
        "Schotten-Baumann to ester",
    ]

    # List of hydrolysis reactions (reverse of esterification)
    hydrolysis_reactions = [
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Ester saponification (methyl deprotection)",
        "Ester saponification (alkyl deprotection)",
        "COOH ethyl deprotection",
    ]

    def dfs_traverse(node):
        nonlocal esterification_count

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an esterification reaction using the checker
            is_esterification = False

            # Check direct esterification reactions
            for reaction_type in esterification_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    is_esterification = True
                    print(f"Detected esterification reaction: {reaction_type}")
                    break

            # Check hydrolysis reactions (which are esterifications in retrosynthesis)
            if not is_esterification:
                for reaction_type in hydrolysis_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_esterification = True
                        print(
                            f"Detected esterification reaction (reverse of hydrolysis): {reaction_type}"
                        )
                        break

            # If not directly identified, check for functional group conversions
            if not is_esterification:
                # Check if reactants have carboxylic acid and product has ester
                reactant_has_acid = False
                reactant_has_alcohol = False
                reactant_has_acyl_halide = False
                reactant_has_ester = False

                for reactant in reactants:
                    if checker.check_fg("Carboxylic acid", reactant):
                        reactant_has_acid = True
                    if any(
                        checker.check_fg(alcohol_type, reactant)
                        for alcohol_type in [
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Aromatic alcohol",
                            "Enol",
                        ]
                    ):
                        reactant_has_alcohol = True
                    if checker.check_fg("Acyl halide", reactant):
                        reactant_has_acyl_halide = True
                    if checker.check_fg("Ester", reactant):
                        reactant_has_ester = True

                product_has_ester = checker.check_fg("Ester", product)
                product_has_acid = checker.check_fg("Carboxylic acid", product)

                # Verify acid is converted to ester (not just both present)
                if product_has_ester and (reactant_has_acid or reactant_has_acyl_halide):
                    # Count acids in reactants and product to ensure conversion
                    acid_count_reactants = sum(
                        1 for r in reactants if checker.check_fg("Carboxylic acid", r)
                    )
                    acid_count_product = 1 if checker.check_fg("Carboxylic acid", product) else 0

                    if acid_count_reactants > acid_count_product:
                        print("Detected esterification reaction by acid to ester conversion")
                        is_esterification = True

                # Check for alcohol to ester conversion (with acyl compounds)
                elif product_has_ester and reactant_has_alcohol:
                    # Check if any reactant has an acylating agent
                    has_acylating_agent = reactant_has_acyl_halide or any(
                        checker.check_fg("Anhydride", r) for r in reactants
                    )

                    if has_acylating_agent:
                        print("Detected esterification reaction by alcohol acylation")
                        is_esterification = True

                # Check for transesterification (ester to different ester)
                elif product_has_ester and reactant_has_ester and reactant_has_alcohol:
                    # Count esters in reactants and product
                    ester_count_reactants = sum(
                        1 for r in reactants if checker.check_fg("Ester", r)
                    )

                    # If number of esters didn't increase, it's likely a transesterification
                    if ester_count_reactants >= 1:
                        print("Detected transesterification (ester to different ester)")
                        is_esterification = True

                # Check for ester hydrolysis (which is esterification in retrosynthesis)
                elif product_has_acid and reactant_has_ester:
                    # Count esters in reactants and product
                    ester_count_reactants = sum(
                        1 for r in reactants if checker.check_fg("Ester", r)
                    )
                    ester_count_product = 1 if checker.check_fg("Ester", product) else 0

                    if ester_count_reactants > ester_count_product:
                        print("Detected esterification (reverse of ester hydrolysis)")
                        is_esterification = True

            if is_esterification:
                esterification_count += 1
                print(f"Current esterification count: {esterification_count}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Total esterification reactions found: {esterification_count}")
    return esterification_count >= 2
