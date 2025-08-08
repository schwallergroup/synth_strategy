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
    This function detects the use of THP (tetrahydropyran) protection strategy in the synthesis.
    It looks for THP protection events throughout the route.
    """
    thp_protection_count = 0

    def dfs_traverse(node):
        nonlocal thp_protection_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if product has THP
                product_has_thp = checker.check_ring("tetrahydropyran", product)

                # Check if any reactant has THP
                reactants_with_thp = sum(
                    1 for reactant in reactants if checker.check_ring("tetrahydropyran", reactant)
                )

                # If product has THP but not all reactants have THP, it's a protection event
                if product_has_thp and reactants_with_thp < len(reactants):
                    # Make sure at least one reactant doesn't have THP
                    for reactant in reactants:
                        if not checker.check_ring("tetrahydropyran", reactant):
                            thp_protection_count += 1
                            print(f"THP protection detected at reaction: {rsmi}")
                            break

                # Also check for alcohol protection specifically
                if product_has_thp:
                    for reactant in reactants:
                        if checker.check_fg("Primary alcohol", reactant) or checker.check_fg(
                            "Secondary alcohol", reactant
                        ):
                            if not checker.check_ring("tetrahydropyran", reactant):
                                # This is an alcohol being protected with THP
                                print(f"Alcohol THP protection detected at reaction: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total THP protection events: {thp_protection_count}")

    # Return True if we found at least one THP protection event
    return thp_protection_count >= 1
