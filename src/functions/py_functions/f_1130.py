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
    This function detects a strategy where a silyl protection group (TBDMS) is installed
    early and maintained throughout most of the synthesis.
    """
    protection_step_found = False
    protection_maintained = False
    protected_mol_at_depth = {}  # Track protected molecules at each depth

    # Check if starting material is already protected
    def check_starting_materials(node):
        if node["type"] == "mol" and node.get("in_stock", False):
            if checker.check_fg(
                "TMS ether protective group", node["smiles"]
            ) or checker.check_fg("Silyl protective group", node["smiles"]):
                return True
        for child in node.get("children", []):
            if check_starting_materials(child):
                return True
        return False

    starting_material_protected = check_starting_materials(route)
    if starting_material_protected:
        print("Starting material already contains silyl protection")

    def dfs_traverse(node, current_depth=0):
        nonlocal protection_step_found, protection_maintained

        # Extract depth from ID if available
        node_depth = current_depth
        if node.get("metadata", {}).get("ID", "") and "Depth: " in node.get(
            "metadata", {}
        ).get("ID", ""):
            try:
                depth_str = (
                    node.get("metadata", {})
                    .get("ID", "")
                    .split("Depth: ")[1]
                    .split()[0]
                )
                node_depth = int(depth_str)
            except (ValueError, IndexError):
                pass

        if node["type"] == "mol":
            # Check if molecule has silyl protection
            mol_smiles = node["smiles"]
            if checker.check_fg(
                "TMS ether protective group", mol_smiles
            ) or checker.check_fg("Silyl protective group", mol_smiles):
                protected_mol_at_depth[node_depth] = mol_smiles
                print(f"Detected protected molecule at depth {node_depth}")

                # If we see protection at a late stage (low depth)
                if node_depth <= 3:
                    protection_maintained = True
                    print(
                        f"Confirmed silyl protection at late stage (depth {node_depth})"
                    )

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains silyl protected oxygen
            product_has_silyl = checker.check_fg(
                "TMS ether protective group", product_smiles
            ) or checker.check_fg("Silyl protective group", product_smiles)

            # Check if any reactant already has silyl protection
            reactants_have_silyl = any(
                checker.check_fg("TMS ether protective group", r)
                or checker.check_fg("Silyl protective group", r)
                for r in reactants_smiles
            )

            # Check if this is a silyl protection reaction
            is_protection_reaction = checker.check_reaction(
                "Alcohol protection with silyl ethers", rsmi
            )

            # If this is a protection step (early in synthesis)
            if is_protection_reaction or (
                product_has_silyl and not reactants_have_silyl
            ):
                protection_step_found = True
                protected_mol_at_depth[node_depth] = product_smiles
                print(f"Detected silyl protection installation at depth {node_depth}")

            # If protection is maintained in subsequent steps
            elif product_has_silyl and reactants_have_silyl:
                protected_mol_at_depth[node_depth] = product_smiles
                print(f"Detected silyl protection maintained at depth {node_depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, node_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Verify we have evidence of protection strategy
    has_protection_strategy = False

    # Case 1: Protection step found and maintained to late stages
    if protection_step_found and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        if len(depths) >= 2:
            has_protection_strategy = True
            print(
                f"Detected silyl protection strategy from depth {max(depths)} to {min(depths)}"
            )

    # Case 2: Starting material was already protected and maintained
    elif starting_material_protected and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        if len(depths) >= 2:
            has_protection_strategy = True
            print(
                f"Detected silyl protection maintained from starting material to late stage"
            )

    # Case 3: Protection is observed at multiple stages including late stage
    elif len(protected_mol_at_depth) >= 2 and protection_maintained:
        depths = sorted(protected_mol_at_depth.keys())
        has_protection_strategy = True
        print(f"Detected silyl protection at multiple stages including late stage")

    return has_protection_strategy
