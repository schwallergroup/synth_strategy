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
    Detects if the synthesis uses Weinreb amide as an intermediate for carbonyl activation,
    followed by nucleophilic addition (e.g., Grignard) to form a ketone.
    """
    # Track Weinreb amides and their conversions
    weinreb_amides = set()  # Store SMILES of found Weinreb amides
    weinreb_to_ketone_reactions = []  # Store reactions converting Weinreb amides to ketones

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            # Check if this molecule is a Weinreb amide
            mol_smiles = node["smiles"]
            if checker.check_fg("Weinreb amide", mol_smiles):
                weinreb_amides.add(mol_smiles)
                print(f"Found Weinreb amide at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant is a Weinreb amide
            weinreb_reactant = None
            for reactant in reactants:
                if checker.check_fg("Weinreb amide", reactant):
                    weinreb_reactant = reactant
                    break

                # Fallback check for tertiary amides with methoxy group
                if (
                    checker.check_fg("Tertiary amide", reactant)
                    and "OC" in reactant
                    and "N(C)" in reactant
                ):
                    weinreb_reactant = reactant
                    weinreb_amides.add(reactant)
                    print(
                        f"Found potential Weinreb amide (tertiary amide with methoxy) at depth {depth}: {reactant}"
                    )
                    break

            # If a Weinreb amide is in the reactants, check if product is a ketone
            if weinreb_reactant and checker.check_fg("Ketone", product):
                print(f"Found conversion of Weinreb amide to ketone at depth {depth}: {rsmi}")
                weinreb_to_ketone_reactions.append(rsmi)

            # Check for specific reaction types
            if checker.check_reaction("Ketone from Weinreb amide", rsmi):
                print(f"Found specific Weinreb amide to ketone reaction at depth {depth}: {rsmi}")
                weinreb_to_ketone_reactions.append(rsmi)

            # Check for Grignard reaction with Weinreb amide
            if any(checker.check_fg("Weinreb amide", r) for r in reactants) and any(
                checker.check_fg("Magnesium halide", r) for r in reactants
            ):
                if checker.check_fg("Ketone", product):
                    print(f"Found Grignard reaction with Weinreb amide at depth {depth}: {rsmi}")
                    weinreb_to_ketone_reactions.append(rsmi)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we found both Weinreb amide and its conversion to ketone
    result = len(weinreb_amides) > 0 and len(weinreb_to_ketone_reactions) > 0
    print(
        f"Weinreb amides found: {len(weinreb_amides)}, Conversions to ketone: {len(weinreb_to_ketone_reactions)}"
    )
    return result
