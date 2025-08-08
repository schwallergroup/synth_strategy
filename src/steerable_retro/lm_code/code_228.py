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
    This function detects preservation of heterocyclic rings (chloropyridine and piperidine)
    throughout the synthesis.
    """
    # Initialize tracking variables
    target_has_chloropyridine = False
    target_has_piperidine = False

    # First check if the target molecule contains the heterocycles
    if route["type"] == "mol":
        target_smiles = route["smiles"]
        target_has_chloropyridine = checker.check_ring(
            "pyridine", target_smiles
        ) and checker.check_fg("Aromatic halide", target_smiles)
        target_has_piperidine = checker.check_ring("piperidine", target_smiles)

        print(f"Target molecule has chloropyridine: {target_has_chloropyridine}")
        print(f"Target molecule has piperidine: {target_has_piperidine}")

        # If neither heterocycle is present in the target, return True (nothing to preserve)
        if not target_has_chloropyridine and not target_has_piperidine:
            print("No heterocycles to preserve in target molecule")
            return True

    # Track preservation through synthesis
    def dfs_traverse(node, depth=0, has_chloropyridine=False, has_piperidine=False):
        # For molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for heterocycles in this molecule
            current_has_chloropyridine = checker.check_ring(
                "pyridine", mol_smiles
            ) and checker.check_fg("Aromatic halide", mol_smiles)
            current_has_piperidine = checker.check_ring("piperidine", mol_smiles)

            # Update tracking for this branch
            branch_has_chloropyridine = has_chloropyridine or current_has_chloropyridine
            branch_has_piperidine = has_piperidine or current_has_piperidine

            # If this is a starting material, check if we've preserved what we need to
            if node.get("in_stock", False):
                # If target has chloropyridine but this branch doesn't, preservation failed
                if target_has_chloropyridine and not branch_has_chloropyridine:
                    print(f"Chloropyridine not preserved in branch ending with {mol_smiles}")
                    return False

                # If target has piperidine but this branch doesn't, preservation failed
                if target_has_piperidine and not branch_has_piperidine:
                    print(f"Piperidine not preserved in branch ending with {mol_smiles}")
                    return False

                # This branch preserves needed heterocycles
                return True

            # For non-stock molecules, continue traversal with updated tracking
            for child in node.get("children", []):
                if not dfs_traverse(
                    child, depth + 1, branch_has_chloropyridine, branch_has_piperidine
                ):
                    return False

            return True

        # For reaction nodes
        elif node["type"] == "reaction":
            # Get reaction SMILES if available
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check if heterocycles are preserved in this reaction
                if has_chloropyridine and not (
                    checker.check_ring("pyridine", reactants)
                    and checker.check_fg("Aromatic halide", reactants)
                ):
                    print(f"Chloropyridine lost in reaction: {rsmi}")
                    return False

                if has_piperidine and not checker.check_ring("piperidine", reactants):
                    print(f"Piperidine lost in reaction: {rsmi}")
                    return False

            # Continue traversal
            for child in node.get("children", []):
                if not dfs_traverse(child, depth + 1, has_chloropyridine, has_piperidine):
                    return False

            return True

        return True

    # Start traversal from the target molecule
    result = dfs_traverse(route, 0, target_has_chloropyridine, target_has_piperidine)

    if result:
        print("Heterocycles preserved throughout synthesis")
    else:
        print("Heterocycles not preserved throughout synthesis")

    return result
