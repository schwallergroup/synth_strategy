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
    Detects a synthesis strategy involving:
    1. BOC protection maintained throughout synthesis
    2. Heterocycle construction/functionalization
    3. No protection/deprotection steps
    """
    # List of common heterocycles to check
    heterocycle_types = [
        "pyrazole",
        "pyrrole",
        "pyridine",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "pyrimidine",
        "piperidine",
        "morpholine",
        "piperazine",
        "indole",
        "benzimidazole",
        "quinoline",
        "isoquinoline",
    ]

    # Track if we've found a heterocycle
    has_heterocycle = False

    # Track if we've found any Boc protection/deprotection reactions
    has_boc_protection_deprotection = False

    # Track synthesis paths and BOC presence
    paths_with_boc = []

    def dfs_traverse(node, current_path=None, depth=0):
        nonlocal has_heterocycle, has_boc_protection_deprotection

        if current_path is None:
            current_path = []

        # Add current node to path
        current_path = current_path + [node]

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for heterocycles
            if not has_heterocycle:
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, mol_smiles):
                        print(f"Found heterocycle: {heterocycle} in molecule: {mol_smiles}")
                        has_heterocycle = True
                        break

            # If this is a leaf node (starting material), check if it has BOC
            if node.get("in_stock", False) or not node.get("children", []):
                has_boc = checker.check_fg("Boc", mol_smiles)
                if has_boc:
                    print(f"Found Boc in starting material: {mol_smiles}")
                    paths_with_boc.append((current_path, True))  # Path starts with BOC
                else:
                    paths_with_boc.append((current_path, False))  # Path starts without BOC

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]

                # Check if this is a Boc protection or deprotection reaction
                if (
                    checker.check_reaction("Boc amine protection", rxn_smiles)
                    or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rxn_smiles)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rxn_smiles)
                    or checker.check_reaction("Boc amine protection of secondary amine", rxn_smiles)
                    or checker.check_reaction("Boc amine protection of primary amine", rxn_smiles)
                    or checker.check_reaction("Boc amine deprotection", rxn_smiles)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rxn_smiles)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rxn_smiles)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rxn_smiles)
                ):
                    print(f"Found Boc protection/deprotection reaction: {rxn_smiles}")
                    has_boc_protection_deprotection = True

        # Process children
        if not node.get("children", []):
            # This is a leaf node in a non-starting material, check if it has BOC
            if node["type"] == "mol" and not node.get("in_stock", False):
                has_boc = checker.check_fg("Boc", node["smiles"])
                paths_with_boc.append((current_path, has_boc))
        else:
            for child in node.get("children", []):
                dfs_traverse(child, current_path, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Process paths to check if BOC is maintained throughout
    all_paths_have_boc = True

    for path, starts_with_boc in paths_with_boc:
        if not starts_with_boc:
            continue  # Skip paths that don't start with BOC

        # Check if BOC is maintained throughout this path
        path_maintains_boc = True
        for node in path:
            if node["type"] == "mol":
                if not checker.check_fg("Boc", node["smiles"]):
                    path_maintains_boc = False
                    break

        if not path_maintains_boc:
            all_paths_have_boc = False
            break

    # Check if the strategy is present
    strategy_present = (
        has_heterocycle and all_paths_have_boc and not has_boc_protection_deprotection
    )

    if strategy_present:
        print("Detected BOC-protected heterocycle strategy")
    else:
        if not has_heterocycle:
            print("No heterocycle found in the synthesis")
        if not all_paths_have_boc:
            print("BOC protection not maintained throughout synthesis")
        if has_boc_protection_deprotection:
            print("Found BOC protection/deprotection steps")

    return strategy_present
