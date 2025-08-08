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
    This function detects a strategy involving isoxazole ring formation
    through the coupling of two fragments, typically involving a chloro-isoxazole
    precursor and an alkyne-containing fragment.
    """
    # Track if we found the strategy
    found_strategy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        # Look at reaction nodes
        if node["type"] == "reaction":
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (convergent step)
                if len(reactants_smiles) >= 2:
                    # Check for isoxazole in product
                    if checker.check_ring("isoxazole", product_smiles):
                        print(f"Found isoxazole in product at depth {depth}")

                        # Check if isoxazole is not present in all reactants
                        isoxazole_in_all_reactants = all(
                            checker.check_ring("isoxazole", smi) for smi in reactants_smiles
                        )

                        if not isoxazole_in_all_reactants:
                            print(f"Isoxazole not present in all reactants at depth {depth}")

                            # Check for alkyne in reactants
                            alkyne_found = False
                            for smi in reactants_smiles:
                                if checker.check_fg("Alkyne", smi):
                                    alkyne_found = True
                                    print(f"Found alkyne-containing fragment at depth {depth}")
                                    break

                            # Check for chloro-containing reactant or other potential isoxazole precursors
                            chloro_or_precursor_found = False
                            for smi in reactants_smiles:
                                # Check for chloro groups
                                if (
                                    checker.check_fg("Primary halide", smi)
                                    or checker.check_fg("Secondary halide", smi)
                                    or checker.check_fg("Tertiary halide", smi)
                                    or checker.check_fg("Aromatic halide", smi)
                                ):
                                    chloro_or_precursor_found = True
                                    print(f"Found halide-containing fragment at depth {depth}")
                                    break
                                # Check for other common isoxazole precursors like oximes, nitriles, etc.
                                if (
                                    checker.check_fg("Nitrile", smi)
                                    or checker.check_fg("Oxime", smi)
                                    or checker.check_fg("Nitro group", smi)
                                    or checker.check_fg("Azide", smi)
                                ):
                                    chloro_or_precursor_found = True
                                    print(f"Found isoxazole precursor fragment at depth {depth}")
                                    break

                            # Check if this is a cycloaddition or related reaction
                            is_cycloaddition = (
                                checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                                or checker.check_reaction(
                                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                                )
                                or checker.check_reaction(
                                    "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                                )
                            )

                            if is_cycloaddition:
                                print(f"Found cycloaddition reaction at depth {depth}")

                            # If we have evidence of convergent isoxazole formation
                            if (alkyne_found and chloro_or_precursor_found) or is_cycloaddition:
                                print(
                                    f"Found convergent isoxazole formation strategy at depth {depth}"
                                )
                                found_strategy = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_strategy
