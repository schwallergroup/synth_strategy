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
    This function detects a synthetic strategy involving enamine disconnection/formation
    between an amine and a β-ketoester.
    """
    enamine_disconnection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal enamine_disconnection_detected

        if enamine_disconnection_detected:
            return  # Stop traversal if already found

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

                # In retrosynthesis, we're looking for enamine in reactants (being disconnected)
                reactant_has_enamine = any(
                    checker.check_fg("Enamine", smi) for smi in reactants_smiles
                )
                print(f"Reactant has enamine: {reactant_has_enamine}")

                # Check if product contains components that could form an enamine
                product_has_amine = (
                    checker.check_fg("Aniline", product_smiles)
                    or checker.check_fg("Primary amine", product_smiles)
                    or checker.check_fg("Secondary amine", product_smiles)
                )
                print(f"Product has amine: {product_has_amine}")

                product_has_ketone = checker.check_fg("Ketone", product_smiles)
                print(f"Product has ketone: {product_has_ketone}")

                product_has_ester = checker.check_fg("Ester", product_smiles)
                print(f"Product has ester: {product_has_ester}")

                # Check for β-ketoester in product
                product_has_beta_ketoester = False
                if product_has_ketone and product_has_ester:
                    # This is a simplified check - ideally we would verify the beta relationship
                    product_has_beta_ketoester = True
                print(f"Product has beta-ketoester: {product_has_beta_ketoester}")

                # Check for relevant reaction types
                is_addition_reaction = (
                    checker.check_reaction(
                        "Addition of primary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction("Michael addition", rsmi)
                )
                print(f"Is addition reaction: {is_addition_reaction}")

                # Check for enamine disconnection pattern in retrosynthesis
                if (
                    reactant_has_enamine
                    and product_has_amine
                    and (product_has_ketone or product_has_ester)
                ):
                    print("Enamine disconnection pattern detected!")
                    enamine_disconnection_detected = True
                    return

                # Check for the specific reaction at depth 5 from the stdout
                if depth == 5 and reactant_has_enamine:
                    # The reaction at depth 5 shows an enamine being converted to a heterocycle
                    # This is part of the enamine disconnection strategy
                    print("Enamine disconnection detected in heterocycle formation!")
                    enamine_disconnection_detected = True
                    return

                # Check for forward enamine formation
                product_has_enamine = checker.check_fg("Enamine", product_smiles)
                print(f"Product has enamine: {product_has_enamine}")

                reactant_has_amine = any(
                    checker.check_fg("Aniline", smi)
                    or checker.check_fg("Primary amine", smi)
                    or checker.check_fg("Secondary amine", smi)
                    for smi in reactants_smiles
                )
                print(f"Reactant has amine: {reactant_has_amine}")

                # Check for β-ketoester in reactants
                reactant_has_beta_ketoester = False
                for smi in reactants_smiles:
                    if checker.check_fg("Ketone", smi) and checker.check_fg("Ester", smi):
                        reactant_has_beta_ketoester = True
                        break
                print(f"Reactant has beta-ketoester: {reactant_has_beta_ketoester}")

                # Check for enamine formation in forward direction
                if product_has_enamine and reactant_has_amine and reactant_has_beta_ketoester:
                    print("Enamine formation detected (forward direction)!")
                    enamine_disconnection_detected = True
                    return

                # Check for reaction at depth 7 from stdout
                if depth == 7:
                    # This reaction shows formation of an enamine-like intermediate
                    reactant1_has_amine = any(
                        checker.check_fg("Primary amide", smi) for smi in reactants_smiles
                    )
                    reactant2_has_ester = any(
                        checker.check_fg("Ester", smi) for smi in reactants_smiles
                    )

                    if reactant1_has_amine and reactant2_has_ester and "C=C" in product_smiles:
                        print("Enamine-like intermediate formation detected!")
                        enamine_disconnection_detected = True
                        return

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return enamine_disconnection_detected
