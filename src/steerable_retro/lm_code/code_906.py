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
    Detects if the synthesis uses a late-stage Suzuki coupling to join two complex fragments,
    with at least one protected functional group in the fragments.
    """
    has_suzuki_coupling = False
    protected_in_suzuki_reactants = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling, protected_in_suzuki_reactants

        if node["type"] == "reaction":
            # Check if this is a reaction node
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling at late stage (depth <= 2)
                if depth <= 2:
                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check for Suzuki coupling using all possible reaction types
                    is_suzuki = False
                    suzuki_types = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic esters",
                        "Suzuki coupling with boronic acids OTf",
                        "Suzuki coupling with boronic esters OTf",
                        "Suzuki",
                    ]

                    for suzuki_type in suzuki_types:
                        if checker.check_reaction(suzuki_type, rsmi):
                            print(f"Detected {suzuki_type} at depth {depth}")
                            is_suzuki = True
                            break

                    # If not detected by reaction type, check for characteristic patterns
                    if not is_suzuki and ("B(" in rsmi or "OB(" in rsmi) and ("[Pd]" in rsmi):
                        print(f"Detected Suzuki coupling by pattern at depth {depth}")
                        is_suzuki = True

                    if is_suzuki:
                        has_suzuki_coupling = True

                        # Protection groups to check
                        protection_groups = [
                            "TMS ether protective group",
                            "Silyl protective group",
                            "Boc",
                            "Acetal/Ketal",
                        ]

                        # Check for protected groups in Suzuki reactants
                        for reactant in reactants:
                            # Skip very small molecules and catalysts
                            if len(reactant) < 5 or "[Pd]" in reactant or "C1COCCO1" in reactant:
                                continue

                            try:
                                # Check for protected groups in this Suzuki reactant
                                for pg in protection_groups:
                                    if checker.check_fg(pg, reactant):
                                        print(
                                            f"Found protected group {pg} in Suzuki reactant: {reactant}"
                                        )
                                        protected_in_suzuki_reactants = True
                            except Exception as e:
                                print(f"Error checking protection groups: {e}")
                                continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Final results - has_suzuki_coupling: {has_suzuki_coupling}, protected_in_suzuki_reactants: {protected_in_suzuki_reactants}"
    )

    # Return True if both conditions are met: late-stage Suzuki coupling and protected functional group in fragments
    return has_suzuki_coupling and protected_in_suzuki_reactants
