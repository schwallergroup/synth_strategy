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
    This function detects if a Suzuki coupling is used in the synthesis.
    It looks for C-C bond formation between aromatic rings with boronic acid and aryl halide reactants.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found

        print(f"Traversing node at depth {depth}: {node.get('type', 'unknown')}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction SMILES: {rsmi}")

            # Check each Suzuki reaction type individually
            suzuki_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "{Suzuki}",
            ]

            for suzuki_type in suzuki_types:
                is_suzuki = checker.check_reaction(suzuki_type, rsmi)
                print(f"Checking {suzuki_type}: {is_suzuki}")
                if is_suzuki:
                    suzuki_coupling_found = True
                    print(f"Found Suzuki coupling: {rsmi}")
                    break

            # If not detected by reaction checkers, manually check for Suzuki pattern
            if not suzuki_coupling_found:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for boronic acid or ester in reactants
                    boronic_present = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                        if r
                    )

                    # Check for aryl halide or triflate in reactants
                    halide_present = any(
                        checker.check_fg("Aromatic halide", r) or checker.check_fg("Triflate", r)
                        for r in reactants
                        if r
                    )

                    print(f"Boronic present: {boronic_present}, Halide present: {halide_present}")

                    # Check for Pd catalyst in reagents
                    reagents = rsmi.split(">")[1].split(".")
                    pd_present = any("[Pd]" in r for r in reagents if r)
                    phosphine_present = any("P(" in r for r in reagents if r)

                    print(f"Pd present: {pd_present}, Phosphine present: {phosphine_present}")

                    # If we have boronic acid/ester, aryl halide, and Pd catalyst, it's likely a Suzuki coupling
                    if boronic_present and halide_present and (pd_present or phosphine_present):
                        suzuki_coupling_found = True
                        print(
                            f"Manually identified Suzuki coupling with boronic compound and halide/triflate"
                        )
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    print(f"Suzuki coupling found: {suzuki_coupling_found}")
    return suzuki_coupling_found
