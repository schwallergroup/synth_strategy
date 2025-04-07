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
    Detects a synthetic strategy involving Suzuki coupling using boronic acid intermediates.
    """
    borylation_found = False
    borylation_depth = -1
    suzuki_coupling_found = False
    suzuki_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal borylation_found, borylation_depth, suzuki_coupling_found, suzuki_coupling_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for borylation (preparation of boronic acids/esters)
                if any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in [
                        "Preparation of boronic acids",
                        "Preparation of boronic acids without boronic ether",
                        "Preparation of boronic acids from trifluoroborates",
                        "Preparation of boronic ethers",
                    ]
                ):
                    borylation_found = True
                    borylation_depth = depth
                    print(f"Found borylation at depth {depth}, rsmi: {rsmi}")

                # Alternative check for borylation: aryl halide to boronic acid/ester
                if not borylation_found:
                    aryl_halide_in_reactants = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    )
                    boronic_acid_in_product = checker.check_fg("Boronic acid", product)
                    boronic_ester_in_product = checker.check_fg("Boronic ester", product)

                    if aryl_halide_in_reactants and (
                        boronic_acid_in_product or boronic_ester_in_product
                    ):
                        borylation_found = True
                        borylation_depth = depth
                        print(f"Found borylation (FG check) at depth {depth}, rsmi: {rsmi}")

                # Check for Suzuki coupling
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                ):
                    suzuki_coupling_found = True
                    suzuki_coupling_depth = depth
                    print(f"Found Suzuki coupling at depth {depth}, rsmi: {rsmi}")

                # Alternative check for Suzuki coupling
                if not suzuki_coupling_found:
                    boronic_acid_in_reactants = any(
                        checker.check_fg("Boronic acid", r) for r in reactants
                    )
                    boronic_ester_in_reactants = any(
                        checker.check_fg("Boronic ester", r) for r in reactants
                    )
                    aryl_halide_in_reactants = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    )
                    triflate_in_reactants = any(checker.check_fg("Triflate", r) for r in reactants)

                    if (boronic_acid_in_reactants or boronic_ester_in_reactants) and (
                        aryl_halide_in_reactants or triflate_in_reactants
                    ):
                        suzuki_coupling_found = True
                        suzuki_coupling_depth = depth
                        print(f"Found Suzuki coupling (FG check) at depth {depth}, rsmi: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # In retrosynthetic analysis, borylation should have higher depth than Suzuki coupling
    # (borylation happens earlier in the forward synthesis, so appears deeper in retrosynthetic tree)
    strategy_present = (
        borylation_found
        and suzuki_coupling_found
        and borylation_depth
        > suzuki_coupling_depth  # In retrosynthesis, borylation appears at greater depth
    )

    if strategy_present:
        print(
            f"Found Suzuki coupling strategy with boronic acid intermediate: borylation at depth {borylation_depth}, coupling at depth {suzuki_coupling_depth}"
        )
    else:
        if borylation_found and suzuki_coupling_found:
            print(
                f"Found both reactions but in wrong order: borylation at depth {borylation_depth}, coupling at depth {suzuki_coupling_depth}"
            )
        elif borylation_found:
            print(f"Found only borylation at depth {borylation_depth}")
        elif suzuki_coupling_found:
            print(f"Found only Suzuki coupling at depth {suzuki_coupling_depth}")
        else:
            print("Neither borylation nor Suzuki coupling found")

    return strategy_present
