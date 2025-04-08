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
    This function detects if the synthetic route contains an ester hydrolysis reaction.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Primary check: Use the reaction checker directly
            if checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                print(f"Ester hydrolysis detected via reaction check: {rsmi}")
                ester_hydrolysis_detected = True
            else:
                # Secondary check: Look for ester saponification reactions
                if checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi):
                    print(f"Ester hydrolysis detected via saponification check: {rsmi}")
                    ester_hydrolysis_detected = True
                else:
                    # Fallback: Check for functional group transformation
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    has_ester = any(checker.check_fg("Ester", reactant) for reactant in reactants)
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)

                    # Additional check: Make sure the product doesn't still have an ester
                    # This helps confirm the ester was actually hydrolyzed
                    product_still_has_ester = checker.check_fg("Ester", product)

                    if has_ester and has_carboxylic_acid and not product_still_has_ester:
                        print(
                            f"Ester hydrolysis detected through functional group analysis: {rsmi}"
                        )
                        ester_hydrolysis_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ester_hydrolysis_detected
