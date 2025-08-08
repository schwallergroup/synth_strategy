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
    Detects a linear synthesis with late-stage Suzuki coupling.
    """
    suzuki_coupling_found = False
    linear_synthesis = True

    def dfs_traverse(node, current_depth=0):
        nonlocal suzuki_coupling_found, linear_synthesis

        # Check if this node has more than 2 children (non-linear synthesis)
        if len(node.get("children", [])) > 2:
            linear_synthesis = False
            print(f"Found non-linear branch with {len(node.get('children', []))} children")

        if node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Examining reaction at depth {current_depth}: {rsmi}")

                    # Check for Suzuki coupling using checker function - more comprehensive check
                    is_suzuki = (
                        checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                        or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                        or checker.check_reaction("{Suzuki}", rsmi)
                    )

                    # Fallback check for Suzuki coupling by examining reactants and products
                    if not is_suzuki:
                        # Check for boronic acid/ester and aryl halide in reactants
                        has_boronic = any(
                            "B(O)" in r or "OB(O)" in r or "B(OH)" in r for r in reactants
                        )
                        has_aryl_halide = any(
                            ("Br" in r or "I" in r or "Cl" in r) and any(c in r for c in ["c", "C"])
                            for r in reactants
                        )
                        # If both key components are present, it's likely a Suzuki coupling
                        if has_boronic and has_aryl_halide:
                            print(
                                f"Detected potential Suzuki coupling based on reactants at depth {current_depth}"
                            )
                            is_suzuki = True

                    if is_suzuki:
                        print(f"Found Suzuki coupling at depth {current_depth}")

                    # Check if this is a late stage reaction (depth 0, 1, or 2)
                    if is_suzuki and current_depth <= 2:
                        suzuki_coupling_found = True
                        print(f"Found late-stage Suzuki coupling at depth {current_depth}")

                    # Check if this is a convergent synthesis (more than 2 reactants)
                    significant_reactants = [r for r in reactants if r.strip()]
                    if len(significant_reactants) > 2:
                        linear_synthesis = False
                        print(f"Found convergent step with {len(significant_reactants)} reactants")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    try:
        dfs_traverse(route)
        print(
            f"Suzuki coupling found: {suzuki_coupling_found}, Linear synthesis: {linear_synthesis}"
        )
        return suzuki_coupling_found and linear_synthesis
    except Exception as e:
        print(f"Error traversing route: {e}")
        return False
