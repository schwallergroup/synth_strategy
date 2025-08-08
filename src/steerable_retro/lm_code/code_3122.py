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
    Detects a synthetic strategy involving a Suzuki coupling in the late stage (low depth)
    of the synthesis, specifically looking for boronic acid/ester and aryl halide coupling.
    """
    suzuki_coupling_found = False
    suzuki_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found, suzuki_depth

        if node["type"] == "reaction":
            # Check if this is a Suzuki coupling reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for Suzuki coupling using the checker function
                is_suzuki = (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                    or checker.check_reaction("{Suzuki}", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                )

                if is_suzuki:
                    print(f"Suzuki coupling detected at depth {depth}")
                    suzuki_coupling_found = True
                    suzuki_depth = min(suzuki_depth, depth)
                else:
                    # Manual check for Suzuki coupling components
                    reactants = rsmi.split(">")[0].split(".")

                    # Check for boronic acid/ester in reactants
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    )

                    # Check for aryl halide or triflate in reactants
                    has_aryl_halide = any(
                        checker.check_fg("Aromatic halide", r) or checker.check_fg("Triflate", r)
                        for r in reactants
                    )

                    print(
                        f"Depth {depth} - Is Suzuki: {is_suzuki}, Has boronic: {has_boronic}, Has aryl halide: {has_aryl_halide}"
                    )

                    # Debug reactants
                    for r in reactants:
                        print(f"  Reactant: {r}")
                        print(f"    Boronic acid: {checker.check_fg('Boronic acid', r)}")
                        print(f"    Boronic ester: {checker.check_fg('Boronic ester', r)}")
                        print(f"    Aromatic halide: {checker.check_fg('Aromatic halide', r)}")
                        print(f"    Triflate: {checker.check_fg('Triflate', r)}")

                    # If we have both boronic acid/ester and aryl halide, it's likely a Suzuki coupling
                    if has_boronic and has_aryl_halide:
                        print(f"Manual Suzuki coupling detection at depth {depth}")
                        suzuki_coupling_found = True
                        suzuki_depth = min(suzuki_depth, depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Suzuki coupling found: {suzuki_coupling_found}, at depth: {suzuki_depth}")
    # Consider it late-stage if it occurs at depth 0, 1, 2, or 3
    return suzuki_coupling_found and suzuki_depth <= 3
