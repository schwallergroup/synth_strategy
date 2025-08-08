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
    Detects if the synthesis route involves a reductive amination step.
    """
    has_reductive_amination = False

    def dfs_traverse(node):
        nonlocal has_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a reductive amination reaction using the checker function
            if (
                checker.check_reaction("reductive amination with aldehyde", rsmi)
                or checker.check_reaction("reductive amination with ketone", rsmi)
                or checker.check_reaction("reductive amination with alcohol", rsmi)
                or checker.check_reaction("{reductive amination}", rsmi)
            ):

                print(f"Found reductive amination reaction: {rsmi}")
                has_reductive_amination = True
                return

            # If direct reaction check fails, try to analyze the reaction components
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carbonyl compounds and amines in reactants
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)

                # Check for amine in product
                has_amine_product = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                # If we have the right reactants and product, it might be reductive amination
                if has_amine_product and (
                    (has_aldehyde or has_ketone) and (has_primary_amine or has_secondary_amine)
                ):
                    print(
                        f"Detected potential reductive amination based on functional groups: {rsmi}"
                    )
                    has_reductive_amination = True
            except Exception as e:
                print(f"Error in reductive amination analysis: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_reductive_amination
