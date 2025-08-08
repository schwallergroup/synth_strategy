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
    This function detects a synthetic strategy where a biaryl system is formed via Suzuki coupling
    in a late stage of the synthesis, where one fragment contains a pyrazole heterocycle and
    the other contains a trifluoromethoxy group.
    """
    # Initialize tracking variables
    suzuki_reaction_node = None
    suzuki_depth = float("inf")
    has_trifluoromethoxy = False
    has_pyrazole = False
    all_reaction_nodes = []

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_reaction_node, suzuki_depth, has_trifluoromethoxy, has_pyrazole

        # Check for trifluoromethoxy group in molecule nodes
        if node.get("type") == "mol" and node.get("smiles"):
            mol_smiles = node.get("smiles")
            if checker.check_fg("Trifluoro group", mol_smiles):
                has_trifluoromethoxy = True
                print(f"Found trifluoromethoxy group in molecule: {mol_smiles}")

            if checker.check_ring("pyrazole", mol_smiles):
                has_pyrazole = True
                print(f"Found pyrazole ring in molecule: {mol_smiles}")

        # Check for reactions
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            all_reaction_nodes.append((node, depth))

            # Check for Suzuki coupling at late stage (depth <= 2)
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
            )

            print(f"Checking reaction at depth {depth}, is Suzuki: {is_suzuki}, RSMI: {rsmi}")

            if depth <= 2 and is_suzuki:
                # Only update if this is a later stage (lower depth) than previously found
                if depth < suzuki_depth:
                    suzuki_reaction_node = node
                    suzuki_depth = depth
                    print(f"Found late-stage Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If no late-stage Suzuki found, check for any Suzuki coupling with the required fragments
    if not suzuki_reaction_node:
        print(
            "No late-stage Suzuki found, checking for any Suzuki coupling with required fragments"
        )
        for node, depth in all_reaction_nodes:
            rsmi = node["metadata"]["rsmi"]
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
            )

            if is_suzuki:
                print(f"Found Suzuki coupling at depth {depth}: {rsmi}")
                suzuki_reaction_node = node
                suzuki_depth = depth
                break

    # If still no Suzuki found, check for biaryl formation manually
    if not suzuki_reaction_node:
        print("No Suzuki coupling found, checking for biaryl formation manually")
        for node, depth in all_reaction_nodes:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check if one reactant has pyrazole and another has trifluoromethoxy
            reactant_with_pyrazole = None
            reactant_with_trifluoromethoxy = None

            for r in reactants:
                if checker.check_ring("pyrazole", r):
                    reactant_with_pyrazole = r
                if checker.check_fg("Trifluoro group", r):
                    reactant_with_trifluoromethoxy = r

            # Check if product has both groups
            product_has_both = checker.check_ring("pyrazole", product_part) and checker.check_fg(
                "Trifluoro group", product_part
            )

            # If both groups are found on different reactants and the product has both
            if (
                reactant_with_pyrazole
                and reactant_with_trifluoromethoxy
                and reactant_with_pyrazole != reactant_with_trifluoromethoxy
                and product_has_both
            ):
                print(f"Found biaryl formation at depth {depth}: {rsmi}")
                suzuki_reaction_node = node
                suzuki_depth = depth
                break

    # Verify that pyrazole and trifluoromethoxy are on different fragments in the Suzuki coupling
    correct_fragments = False
    if suzuki_reaction_node:
        rsmi = suzuki_reaction_node["metadata"]["rsmi"]
        reactants_part = rsmi.split(">")[0]
        product_part = rsmi.split(">")[-1]
        reactants = reactants_part.split(".")

        print(f"Analyzing coupling reaction: {rsmi}")
        print(f"Reactants: {reactants}")
        print(f"Product: {product_part}")

        # Check if one reactant has pyrazole and another has trifluoromethoxy
        reactant_with_pyrazole = None
        reactant_with_trifluoromethoxy = None

        for r in reactants:
            if checker.check_ring("pyrazole", r):
                reactant_with_pyrazole = r
                print(f"Found reactant with pyrazole: {r}")
            if checker.check_fg("Trifluoro group", r):
                reactant_with_trifluoromethoxy = r
                print(f"Found reactant with trifluoromethoxy: {r}")

        # Verify that the product contains both groups (biaryl formation)
        product_has_both = checker.check_ring("pyrazole", product_part) and checker.check_fg(
            "Trifluoro group", product_part
        )

        # If both groups are found on different reactants and the product has both
        if (
            reactant_with_pyrazole
            and reactant_with_trifluoromethoxy
            and reactant_with_pyrazole != reactant_with_trifluoromethoxy
            and product_has_both
        ):
            correct_fragments = True
            print(
                "Pyrazole and trifluoromethoxy are on different coupling fragments and joined in product"
            )

    # Check if the strategy is present
    # For late-stage, we'll consider depth <= 3 instead of just 2
    strategy_present = (
        suzuki_depth <= 3 and has_pyrazole and has_trifluoromethoxy and correct_fragments
    )

    if strategy_present:
        print(
            f"Detected biaryl Suzuki coupling with pyrazole and trifluoromethoxy groups at depth {suzuki_depth}"
        )
    else:
        print("Strategy not detected")
        if suzuki_depth > 3:
            print(f"- No late-stage coupling found (depth: {suzuki_depth})")
        elif suzuki_depth == float("inf"):
            print("- No suitable coupling reaction found at all")
        if not has_pyrazole:
            print("- No pyrazole ring found")
        if not has_trifluoromethoxy:
            print("- No trifluoromethoxy group found")
        if not correct_fragments:
            print(
                "- Pyrazole and trifluoromethoxy not on different coupling fragments or not joined in product"
            )

    return strategy_present
