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
    This function detects if the final step involves amide bond formation.
    """
    # Find the final reaction step (closest to the target molecule)
    final_reaction = None

    def find_final_reaction(node, depth=0):
        nonlocal final_reaction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # If this is the first reaction we encounter (lowest depth), it's the final step
            if final_reaction is None or depth < final_reaction["depth"]:
                final_reaction = {"node": node, "depth": depth}

        # Recursively check children
        for child in node.get("children", []):
            find_final_reaction(child, depth + 1)

    # Start traversal from the root
    find_final_reaction(route)

    # If no reaction was found, return False
    if final_reaction is None:
        print("No reaction found in the route")
        return False

    # Extract the final reaction
    final_step = final_reaction["node"]
    rsmi = final_step["metadata"]["rsmi"]
    reactants = rsmi.split(">")[0].split(".")
    product = rsmi.split(">")[-1]

    print(f"Analyzing final step: {rsmi}")

    # Check if the reaction is an amide formation using the checker functions
    amide_formation_reactions = [
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Carboxylic acid with primary amine to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Acyl chloride with ammonia to amide",
        "Ester with ammonia to amide",
        "Schotten-Baumann_amide",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of primary amines",
        "Acylation of secondary amines",
    ]

    # Check if the reaction is any of the amide formation types
    for reaction_type in amide_formation_reactions:
        if checker.check_reaction(reaction_type, rsmi):
            print(f"Detected amide formation reaction: {reaction_type}")
            return True

    # If reaction type check failed, check for functional group changes
    product_mol = Chem.MolFromSmiles(product)
    if not product_mol:
        print("Could not parse product molecule")
        return False

    # Check if product contains an amide group
    has_amide = (
        checker.check_fg("Primary amide", product)
        or checker.check_fg("Secondary amide", product)
        or checker.check_fg("Tertiary amide", product)
    )

    if not has_amide:
        print("Product does not contain an amide group")
        return False

    # Check if reactants contain carboxylic acid and amine
    has_acid = False
    has_amine = False

    for reactant in reactants:
        if checker.check_fg("Carboxylic acid", reactant):
            has_acid = True
            print(f"Found carboxylic acid in reactant: {reactant}")

        if (
            checker.check_fg("Primary amine", reactant)
            or checker.check_fg("Secondary amine", reactant)
            or checker.check_fg("Aniline", reactant)
        ):
            has_amine = True
            print(f"Found amine in reactant: {reactant}")

    # Also check for acyl halides or esters which can form amides
    for reactant in reactants:
        if checker.check_fg("Acyl halide", reactant) or checker.check_fg("Ester", reactant):
            has_acid = True
            print(f"Found acyl halide or ester in reactant: {reactant}")

    # Check if the product has a new amide bond that wasn't in the reactants
    new_amide_formed = False
    for reactant in reactants:
        if (
            checker.check_fg("Primary amide", reactant)
            or checker.check_fg("Secondary amide", reactant)
            or checker.check_fg("Tertiary amide", reactant)
        ):
            # If reactant already has amide, make sure product has more amides
            reactant_mol = Chem.MolFromSmiles(reactant)
            if reactant_mol:
                # Count amide groups in reactant and product
                amide_count_reactant = sum(
                    [
                        len(checker.get_fg_atom_indices("Primary amide", reactant)),
                        len(checker.get_fg_atom_indices("Secondary amide", reactant)),
                        len(checker.get_fg_atom_indices("Tertiary amide", reactant)),
                    ]
                )
                amide_count_product = sum(
                    [
                        len(checker.get_fg_atom_indices("Primary amide", product)),
                        len(checker.get_fg_atom_indices("Secondary amide", product)),
                        len(checker.get_fg_atom_indices("Tertiary amide", product)),
                    ]
                )
                if amide_count_product > amide_count_reactant:
                    new_amide_formed = True
                    print(
                        f"New amide bond formed: {amide_count_reactant} in reactant, {amide_count_product} in product"
                    )

    # If all reactants already had amides and no new ones formed, set new_amide_formed based on acid+amine
    if not new_amide_formed and has_acid and has_amine:
        new_amide_formed = True

    if new_amide_formed:
        print("Amide formation detected based on functional group analysis")
        return True

    print("Not an amide formation reaction")
    return False
