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
    Detects a strategy involving a nitrile intermediate formed from an amine.
    """
    has_nitrile = False
    has_amine_to_nitrile = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile, has_amine_to_nitrile

        if node["type"] == "mol":
            # Check if molecule contains nitrile
            if checker.check_fg("Nitrile", node["smiles"]):
                has_nitrile = True
                print(f"Found nitrile in molecule at depth {depth}: {node['smiles']}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrile in product
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)
                if product_has_nitrile:
                    has_nitrile = True
                    print(f"Found nitrile in product at depth {depth}: {product_smiles}")

                # Check for amine in reactants
                reactant_has_amine = False
                for reactant in reactants_smiles:
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                    ):
                        reactant_has_amine = True
                        print(f"Found amine in reactant at depth {depth}: {reactant}")
                        break

                # Check for nitrile in reactants
                reactant_has_nitrile = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Nitrile", reactant):
                        reactant_has_nitrile = True
                        print(f"Found nitrile in reactant at depth {depth}: {reactant}")
                        break

                # Check for amide in reactants (potential nitrile precursor)
                reactant_has_amide = False
                for reactant in reactants_smiles:
                    if checker.check_fg("Primary amide", reactant) or checker.check_fg(
                        "Secondary amide", reactant
                    ):
                        reactant_has_amide = True
                        print(f"Found amide in reactant at depth {depth}: {reactant}")
                        break

                # In retrosynthesis, we're looking for reactions where nitrile is in the product
                # and amine or amide is in the reactants (which means amine/amide to nitrile in forward direction)
                if (
                    product_has_nitrile
                    and (reactant_has_amine or reactant_has_amide)
                    and not reactant_has_nitrile
                ):
                    print(f"Potential amine/amide to nitrile transformation at depth {depth}")
                    has_amine_to_nitrile = True

                # Check for direct amine to nitrile conversion
                if reactant_has_amine and product_has_nitrile and not reactant_has_nitrile:
                    print(f"Direct amine to nitrile transformation detected at depth {depth}")
                    has_amine_to_nitrile = True

                # Check for nitrile-related reactions
                if product_has_nitrile or reactant_has_nitrile:
                    # Check common reactions involving nitriles
                    if any(
                        [
                            checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi),
                            checker.check_reaction("Reduction of nitrile to amine", rsmi),
                            checker.check_reaction("Nitrile to amide", rsmi),
                            checker.check_reaction("Grignard from nitrile to ketone", rsmi),
                        ]
                    ):
                        print(f"Found nitrile-related reaction at depth {depth}")
                        has_nitrile = True

                # Check for alkene-nitrile addition reactions
                if product_has_nitrile:
                    for reactant in reactants_smiles:
                        if checker.check_fg("Alkene", reactant):
                            print(f"Potential alkene to nitrile transformation at depth {depth}")
                            has_amine_to_nitrile = True
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check for Michael addition with nitrile acceptors
    def check_michael_nitrile(node):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                if (
                    checker.check_reaction("Michael addition", rsmi)
                    or checker.check_reaction("aza-Michael addition primary", rsmi)
                    or checker.check_reaction("aza-Michael addition secondary", rsmi)
                    or checker.check_reaction("aza-Michael addition aromatic", rsmi)
                    or checker.check_reaction("thia-Michael addition", rsmi)
                    or checker.check_reaction("oxa-Michael addition", rsmi)
                ):

                    reactants_smiles = rsmi.split(">")[0].split(".")
                    for reactant in reactants_smiles:
                        if checker.check_fg("Nitrile", reactant) and checker.check_fg(
                            "Alkene", reactant
                        ):
                            print(f"Found Michael addition with nitrile acceptor: {rsmi}")
                            return True

                    # Also check for acrylonitrile which is a common Michael acceptor
                    for reactant in reactants_smiles:
                        if "C=CC#N" in reactant:
                            print(f"Found Michael addition with acrylonitrile: {rsmi}")
                            return True
            except Exception as e:
                print(f"Error checking Michael addition: {e}")

        # Check children
        for child in node.get("children", []):
            if check_michael_nitrile(child):
                return True

        return False

    # Additional check for Michael addition with nitrile acceptors
    michael_nitrile = check_michael_nitrile(route)

    # Check for acrylonitrile in the route
    def check_acrylonitrile(node):
        if node["type"] == "mol":
            if "C=CC#N" in node["smiles"]:
                print(f"Found acrylonitrile: {node['smiles']}")
                return True

        # Check children
        for child in node.get("children", []):
            if check_acrylonitrile(child):
                return True

        return False

    # If we find acrylonitrile, it's a strong indicator of nitrile intermediate strategy
    acrylonitrile_present = check_acrylonitrile(route)
    if acrylonitrile_present:
        has_amine_to_nitrile = True

    strategy_present = has_nitrile and (
        has_amine_to_nitrile or michael_nitrile or acrylonitrile_present
    )

    if strategy_present:
        print("Detected nitrile intermediate strategy")
    else:
        if has_nitrile:
            print("Found nitrile but no amine-to-nitrile transformation")
        else:
            print("No nitrile found in the route")

    return strategy_present
