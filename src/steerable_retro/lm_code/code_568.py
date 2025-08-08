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
    Detects a synthetic strategy involving triazole formation via isothiocyanate-thioester route.
    The strategy involves:
    1. Starting with an aniline derivative
    2. Converting to isothiocyanate
    3. Forming a thioester intermediate
    4. Constructing a triazole ring using hydrazine
    """
    # Initialize flags to track key intermediates and transformations
    has_aniline = False
    has_isothiocyanate = False
    has_thioester = False
    has_triazole = False

    # Track reaction transformations
    aniline_to_isothiocyanate = False
    isothiocyanate_to_thioester = False
    thioester_to_triazole = False
    hydrazine_involved = False

    # Track molecules for sequence verification
    molecules = {}
    reactions = {}

    def dfs_traverse(node, depth=0):
        nonlocal has_aniline, has_isothiocyanate, has_thioester, has_triazole
        nonlocal aniline_to_isothiocyanate, isothiocyanate_to_thioester, thioester_to_triazole, hydrazine_involved

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Store molecule with its depth for sequence analysis
            molecules[mol_smiles] = depth

            # Check for structural features in molecules
            if checker.check_fg("Aniline", mol_smiles) or (
                checker.check_fg("Primary amine", mol_smiles) and "c" in mol_smiles
            ):
                has_aniline = True
                print(f"Found aniline/aromatic amine at depth {depth}: {mol_smiles}")

            if checker.check_fg("Isothiocyanate", mol_smiles):
                has_isothiocyanate = True
                print(f"Found isothiocyanate at depth {depth}: {mol_smiles}")

            # Check for thioester using proper functional group checkers
            if checker.check_fg("Carbo-thioester", mol_smiles) or checker.check_fg(
                "Thiocarbonyl", mol_smiles
            ):
                has_thioester = True
                print(f"Found thioester intermediate at depth {depth}: {mol_smiles}")

            # Check for triazole ring
            if checker.check_ring("triazole", mol_smiles):
                has_triazole = True
                print(f"Found triazole structure at depth {depth}: {mol_smiles}")

            # Check for hydrazine or hydrazine derivatives
            if "NN" in mol_smiles or checker.check_fg("Hydrazine", mol_smiles):
                hydrazine_involved = True
                print(f"Found hydrazine or derivative at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Analyze reaction transformations
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Store reaction with its depth
            reactions[rsmi] = depth

            # Check for aniline to isothiocyanate transformation
            if any(
                checker.check_fg("Aniline", reactant)
                or (checker.check_fg("Primary amine", reactant) and "c" in reactant)
                for reactant in reactants
            ) and checker.check_fg("Isothiocyanate", product):
                aniline_to_isothiocyanate = True
                print(f"Detected aniline to isothiocyanate transformation at depth {depth}")

            # Check for isothiocyanate to thioester transformation
            if any(checker.check_fg("Isothiocyanate", reactant) for reactant in reactants) and (
                checker.check_fg("Carbo-thioester", product)
                or checker.check_fg("Thiocarbonyl", product)
            ):
                isothiocyanate_to_thioester = True
                print(f"Detected isothiocyanate to thioester transformation at depth {depth}")

            # Check for thioester to triazole transformation
            if any(
                checker.check_fg("Carbo-thioester", reactant)
                or checker.check_fg("Thiocarbonyl", reactant)
                for reactant in reactants
            ) and checker.check_ring("triazole", product):
                thioester_to_triazole = True
                print(f"Detected thioester to triazole transformation at depth {depth}")

            # Check for hydrazine involvement in triazole formation
            if any(
                "NN" in reactant or checker.check_fg("Hydrazine", reactant)
                for reactant in reactants
            ) and checker.check_ring("triazole", product):
                thioester_to_triazole = True
                hydrazine_involved = True
                print(f"Detected hydrazine-mediated triazole formation at depth {depth}")

            # Check for direct isothiocyanate to triazole transformation (might skip thioester intermediate)
            if any(
                checker.check_fg("Isothiocyanate", reactant) for reactant in reactants
            ) and checker.check_ring("triazole", product):
                print(f"Detected direct isothiocyanate to triazole transformation at depth {depth}")
                isothiocyanate_to_thioester = True  # Consider this step satisfied
                thioester_to_triazole = True  # Consider this step satisfied

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the key elements of the strategy are present
    # We need isothiocyanate and triazole at minimum, plus at least one key transformation
    key_intermediates = has_isothiocyanate and (has_triazole or has_thioester)
    key_transformations = (
        aniline_to_isothiocyanate or isothiocyanate_to_thioester or thioester_to_triazole
    )

    # For a complete strategy, we'd want to see more elements
    complete_strategy = (
        has_aniline
        and has_isothiocyanate
        and has_triazole
        and (aniline_to_isothiocyanate or isothiocyanate_to_thioester or thioester_to_triazole)
    )

    # Allow for partial strategy detection
    strategy_present = key_intermediates and key_transformations

    if complete_strategy:
        print("Detected complete triazole via isothiocyanate-thioester strategy")
    elif strategy_present:
        print("Detected partial triazole via isothiocyanate strategy")
    else:
        print(
            "Strategy not fully detected. Found: "
            + f"aniline={has_aniline}, isothiocyanate={has_isothiocyanate}, "
            + f"thioester={has_thioester}, triazole={has_triazole}, "
            + f"hydrazine={hydrazine_involved}, "
            + f"aniline→isothiocyanate={aniline_to_isothiocyanate}, "
            + f"isothiocyanate→thioester={isothiocyanate_to_thioester}, "
            + f"thioester→triazole={thioester_to_triazole}"
        )

    return strategy_present or complete_strategy
