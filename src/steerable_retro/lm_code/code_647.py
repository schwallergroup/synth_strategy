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
    This function detects reductive amination as a key fragment coupling strategy.
    It looks for C-N bond formation between an aldehyde-containing fragment and an amine.
    """
    reductive_amination_found = False

    def dfs_traverse(node):
        nonlocal reductive_amination_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction: {rsmi}")

                # Check for reductive amination pattern using checker functions
                if checker.check_reaction("reductive amination with aldehyde", rsmi):
                    print("Reductive amination with aldehyde detected")
                    reductive_amination_found = True
                elif checker.check_reaction("reductive amination with ketone", rsmi):
                    print("Reductive amination with ketone detected")
                    reductive_amination_found = True
                elif checker.check_reaction("reductive amination with alcohol", rsmi):
                    print("Reductive amination with alcohol detected")
                    reductive_amination_found = True
                else:
                    # Manual check for reductive amination pattern
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for aldehyde in reactants
                        aldehyde_present = any(
                            checker.check_fg("Aldehyde", reactant) for reactant in reactants
                        )

                        # Check for primary or secondary amine in reactants
                        amine_present = any(
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            for reactant in reactants
                        )

                        # Check for C-N bond formation (tertiary amine in product)
                        tertiary_amine_in_product = checker.check_fg(
                            "Tertiary amine", product
                        ) or checker.check_fg("Secondary amine", product)

                        if aldehyde_present and amine_present and tertiary_amine_in_product:
                            print(
                                "Manual detection: Reductive amination pattern found (aldehyde + amine)"
                            )
                            reductive_amination_found = True

                        # Also check for ketone-based reductive amination
                        ketone_present = any(
                            checker.check_fg("Ketone", reactant) for reactant in reactants
                        )
                        if ketone_present and amine_present and tertiary_amine_in_product:
                            print(
                                "Manual detection: Reductive amination pattern found (ketone + amine)"
                            )
                            reductive_amination_found = True

                        # Check for reductive amination pattern in reaction name
                        if "{reductive amination}" in rsmi.lower():
                            print("Reductive amination detected in reaction name")
                            reductive_amination_found = True
                    except Exception as e:
                        print(f"Error in manual reductive amination check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {reductive_amination_found}")

    return reductive_amination_found
