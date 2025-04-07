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
    Detects a linear synthesis on a pyridine scaffold with late-stage thioether formation.
    """
    # Track if we found a pyridine scaffold
    has_pyridine = False
    # Track if we found a late-stage thioether formation
    has_late_stage_thioether = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyridine, has_late_stage_thioether

        if node["type"] == "mol":
            # Check if molecule contains pyridine
            if node.get("smiles"):
                # Use checker function to detect pyridine ring
                if checker.check_ring("pyridine", node["smiles"]):
                    has_pyridine = True
                    print(f"Found pyridine scaffold in molecule: {node['smiles']}")

        elif node["type"] == "reaction" and depth <= 2:  # Late stage (depth 0-2)
            # Check for thioether formation in late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if product contains thioether (monosulfide) and pyridine
                if checker.check_fg("Monosulfide", product) and checker.check_ring(
                    "pyridine", product
                ):
                    print(f"Product contains thioether and pyridine: {product}")

                    # Check if any reactant contains thiol
                    has_thiol_reactant = any(
                        checker.check_fg("Aliphatic thiol", reactant)
                        or checker.check_fg("Aromatic thiol", reactant)
                        for reactant in reactants
                    )

                    # Check if this is a thioether formation reaction
                    is_thioether_formation = (
                        checker.check_reaction("S-alkylation of thiols", rsmi)
                        or checker.check_reaction(
                            "S-alkylation of thiols (ethyl)", rsmi
                        )
                        or checker.check_reaction(
                            "S-alkylation of thiols with alcohols", rsmi
                        )
                        or checker.check_reaction(
                            "S-alkylation of thiols with alcohols (ethyl)", rsmi
                        )
                        or checker.check_reaction("thia-Michael addition", rsmi)
                    )

                    if is_thioether_formation:
                        print(f"Confirmed thioether formation reaction: {rsmi}")

                    # Check if thioether is newly formed (not present in all reactants)
                    all_reactants_have_thioether = all(
                        checker.check_fg("Monosulfide", reactant)
                        for reactant in reactants
                    )

                    # If we have a thiol reactant, a thioether formation reaction type, and not all reactants
                    # already have thioether, then this is a late-stage thioether formation
                    if (
                        has_thiol_reactant
                        and is_thioether_formation
                        and not all_reactants_have_thioether
                    ):
                        has_late_stage_thioether = True
                        print(
                            f"Confirmed late-stage thioether formation at depth {depth}"
                        )

                    # Alternative detection: if product has thioether but not all reactants do
                    elif not all_reactants_have_thioether:
                        # Check if any reactant has thiol
                        if has_thiol_reactant:
                            print(
                                f"Detected thioether formation by thiol reactant at depth {depth}"
                            )
                            has_late_stage_thioether = True
                        # Check for nucleophilic substitution with thiol
                        elif any(
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Aromatic halide", reactant)
                            for reactant in reactants
                        ):
                            print(
                                f"Detected thioether formation by nucleophilic substitution at depth {depth}"
                            )
                            has_late_stage_thioether = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Pyridine scaffold: {has_pyridine}, Late-stage thioether: {has_late_stage_thioether}"
    )
    return has_pyridine and has_late_stage_thioether
