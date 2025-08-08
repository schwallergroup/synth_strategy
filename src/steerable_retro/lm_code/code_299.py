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
    This function detects a synthetic strategy involving elaboration of a heterocyclic core
    (pyridine) through sequential transformations without modifying the core structure.
    """
    # Track if we have a heterocyclic core elaboration strategy
    has_heterocyclic_core = False
    elaboration_reactions = 0
    total_reactions = 0

    # Define heterocyclic cores to check
    heterocycles = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazine",
        "imidazole",
        "pyrazole",
        "oxazole",
        "thiazole",
        "triazole",
        "furan",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # Find the target molecule's heterocyclic core
    target_core = None
    if route["type"] == "mol" and "smiles" in route:
        for core in heterocycles:
            if checker.check_ring(core, route["smiles"]):
                target_core = core
                print(f"Target molecule contains {core} core: {route['smiles']}")
                break

    if target_core is None:
        print("Target molecule does not contain a heterocyclic core")
        return False

    def is_elaboration_reaction(rsmi, core_type):
        """Check if a reaction preserves the heterocyclic core while modifying substituents"""
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # In retrosynthetic context, the product is what we're making
        # and reactants are the precursors

        # Check if product has the specific heterocyclic core
        if not checker.check_ring(core_type, product):
            print(f"Product does not contain {core_type} core: {product}")
            return False

        # Check if at least one reactant has the same heterocyclic core
        reactant_has_same_core = False
        for reactant in reactants:
            if checker.check_ring(core_type, reactant):
                reactant_has_same_core = True
                break

        if not reactant_has_same_core:
            print(f"No reactant contains {core_type} core")
            return False

        # If we reach here, both product and at least one reactant have the same heterocyclic core
        # This is the basic requirement for an elaboration reaction

        # Check if this is an elaboration reaction (not modifying the core)
        elaboration_reaction_types = [
            "Suzuki",
            "Buchwald-Hartwig",
            "N-arylation",
            "Sonogashira",
            "Heck",
            "Negishi",
            "Stille",
            "Acylation",
            "Alkylation",
            "Reductive amination",
            "Mitsunobu",
            "Chan-Lam",
            "Methylation",
            "Amidation",
            "Esterification",
            "Etherification",
            "Williamson ether",
            "Schotten-Baumann",
            "Friedel-Crafts",
            "Reductive amination with aldehyde",
            "Reductive amination with ketone",
            "Reductive amination with alcohol",
        ]

        is_elaboration = False
        for rxn_type in elaboration_reaction_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Found {rxn_type} reaction: {rsmi}")
                is_elaboration = True
                break

        # If we didn't find a specific reaction type, check if the reaction
        # involves common functional group transformations that don't affect the heterocyclic core
        if not is_elaboration:
            # Expanded list of functional group transformations
            fg_transformations = [
                ("Primary halide", "Primary alcohol"),
                ("Primary halide", "Primary amine"),
                ("Secondary halide", "Secondary alcohol"),
                ("Secondary halide", "Secondary amine"),
                ("Tertiary halide", "Tertiary alcohol"),
                ("Tertiary halide", "Tertiary amine"),
                ("Aldehyde", "Primary alcohol"),
                ("Aldehyde", "Ester"),
                ("Aldehyde", "Carboxylic acid"),
                ("Ketone", "Secondary alcohol"),
                ("Carboxylic acid", "Ester"),
                ("Nitrile", "Primary amide"),
                ("Nitrile", "Carboxylic acid"),
                ("Nitro group", "Primary amine"),
                ("Aromatic halide", "Ether"),
                ("Aromatic halide", "Amine"),
                ("Aromatic halide", "Nitrile"),
                ("Aromatic halide", "Ester"),
                ("Aromatic halide", "Carboxylic acid"),
                ("Ester", "Carboxylic acid"),
                ("Ester", "Alcohol"),
                ("Ester", "Amide"),
            ]

            for fg1, fg2 in fg_transformations:
                # Check if product has fg1 and reactant has fg2 or vice versa
                product_has_fg1 = checker.check_fg(fg1, product)
                product_has_fg2 = checker.check_fg(fg2, product)

                for reactant in reactants:
                    reactant_has_fg1 = checker.check_fg(fg1, reactant)
                    reactant_has_fg2 = checker.check_fg(fg2, reactant)

                    if (product_has_fg1 and reactant_has_fg2) or (
                        product_has_fg2 and reactant_has_fg1
                    ):
                        print(f"Found functional group transformation: {fg1} <-> {fg2}")
                        is_elaboration = True
                        break

                if is_elaboration:
                    break

        # If we still haven't identified it as an elaboration, check if the core atoms are preserved
        # Since we've already verified the same core type exists in both reactant and product,
        # this is likely an elaboration reaction even if we can't identify the specific transformation
        if not is_elaboration:
            # Default to True if core is preserved in both reactant and product
            is_elaboration = True
            print(f"Core structure preserved, considering as elaboration: {rsmi}")

        return is_elaboration

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocyclic_core, elaboration_reactions, total_reactions, target_core

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]

            if is_elaboration_reaction(rsmi, target_core):
                elaboration_reactions += 1
                print(
                    f"Found elaboration reaction ({elaboration_reactions}/{total_reactions}): {rsmi}"
                )
            else:
                print(f"Non-elaboration reaction: {rsmi}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Determine if this is a heterocyclic core elaboration strategy
    # Criteria: At least 2 elaboration reactions and majority of reactions are elaborations
    if total_reactions > 0 and elaboration_reactions >= 2:
        elaboration_ratio = elaboration_reactions / total_reactions
        has_heterocyclic_core = elaboration_ratio >= 0.5

    if has_heterocyclic_core:
        print(
            f"Detected heterocyclic core elaboration strategy: {elaboration_reactions}/{total_reactions} reactions"
        )
        return True

    print(
        f"Not a heterocyclic core elaboration strategy: {elaboration_reactions}/{total_reactions} reactions"
    )
    return False
