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
    This function detects a strategy involving O-alkylation reactions,
    particularly focusing on pyridine hydroxyl modifications.
    """
    # Track if we found the pattern
    o_alkylation_count = 0

    def dfs_traverse(node):
        nonlocal o_alkylation_count

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for various O-alkylation reactions
                is_o_alkylation = (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("O-alkylation of amides with diazo compounds", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Alcohol to ether", rsmi)
                    or checker.check_reaction("{Williamson ether}", rsmi)
                )

                # Look for hydroxyl pyridine in reactants
                hydroxyl_pyridine_found = False
                for reactant in reactants:
                    if checker.check_ring("pyridine", reactant) and (
                        checker.check_fg("Phenol", reactant)
                        or checker.check_fg("Aromatic alcohol", reactant)
                        or checker.check_fg("Primary alcohol", reactant)
                        or checker.check_fg("Secondary alcohol", reactant)
                        or checker.check_fg("Tertiary alcohol", reactant)
                    ):
                        hydroxyl_pyridine_found = True
                        print(f"Found hydroxyl pyridine in reactant: {reactant}")
                        break

                # Check for conversion to ether
                if hydroxyl_pyridine_found and checker.check_ring("pyridine", product):
                    # Verify the hydroxyl group was converted to an ether
                    if checker.check_fg("Ether", product) and not (
                        checker.check_fg("Phenol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                        or checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    ):
                        print(f"Confirmed O-alkylation on pyridine: {rsmi}")
                        o_alkylation_count += 1

                # Also check for general O-alkylation patterns even if not a specific named reaction
                elif not is_o_alkylation:
                    for reactant in reactants:
                        if checker.check_ring("pyridine", reactant) and (
                            checker.check_fg("Phenol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            if (
                                checker.check_ring("pyridine", product)
                                and checker.check_fg("Ether", product)
                                and not (
                                    checker.check_fg("Phenol", product)
                                    or checker.check_fg("Aromatic alcohol", product)
                                    or checker.check_fg("Primary alcohol", product)
                                    or checker.check_fg("Secondary alcohol", product)
                                    or checker.check_fg("Tertiary alcohol", product)
                                )
                            ):
                                print(f"Found general O-alkylation pattern on pyridine: {rsmi}")
                                o_alkylation_count += 1
                                break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total O-alkylation reactions found: {o_alkylation_count}")
    # Return True if we found at least one O-alkylation
    return o_alkylation_count >= 1
