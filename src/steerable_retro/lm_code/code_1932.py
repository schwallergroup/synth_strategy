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
    This function detects a reductive amination pattern (C=O → C-NR).
    """
    reductive_amination_found = False

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if the reaction is a reductive amination using the checker function
            if (
                checker.check_reaction("reductive amination with aldehyde", rsmi)
                or checker.check_reaction("reductive amination with ketone", rsmi)
                or checker.check_reaction("reductive amination with alcohol", rsmi)
            ):
                print(f"Reductive amination reaction detected: {rsmi}")
                reductive_amination_found = True
                return

            # If direct reaction check fails, try to infer from functional group changes
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carbonyl groups (aldehyde or ketone) in reactants
            reactant_has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
            reactant_has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
            reactant_has_carbonyl = reactant_has_aldehyde or reactant_has_ketone

            # Check for amines in reactants (primary or secondary)
            reactant_has_primary_amine = any(
                checker.check_fg("Primary amine", r) for r in reactants
            )
            reactant_has_secondary_amine = any(
                checker.check_fg("Secondary amine", r) for r in reactants
            )
            reactant_has_amine = reactant_has_primary_amine or reactant_has_secondary_amine

            # Check for amine in product (could be secondary or tertiary depending on starting amine)
            product_has_secondary_amine = checker.check_fg("Secondary amine", product)
            product_has_tertiary_amine = checker.check_fg("Tertiary amine", product)
            product_has_amine = product_has_secondary_amine or product_has_tertiary_amine

            if reactant_has_carbonyl and reactant_has_amine and product_has_amine:
                print(f"Inferred reductive amination pattern from functional groups: {rsmi}")
                print(f"Reactant has aldehyde: {reactant_has_aldehyde}")
                print(f"Reactant has ketone: {reactant_has_ketone}")
                print(f"Reactant has primary amine: {reactant_has_primary_amine}")
                print(f"Reactant has secondary amine: {reactant_has_secondary_amine}")
                print(f"Product has secondary amine: {product_has_secondary_amine}")
                print(f"Product has tertiary amine: {product_has_tertiary_amine}")
                reductive_amination_found = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return reductive_amination_found
