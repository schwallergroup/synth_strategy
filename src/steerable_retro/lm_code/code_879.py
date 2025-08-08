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
    This function detects synthesis strategies involving multiple ether linkages
    between fragments.
    """
    ether_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal ether_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ether formation reactions directly
                is_ether_formation = (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Chan-Lam etherification", rsmi)
                    or checker.check_reaction("{Williamson ether}", rsmi)
                )

                # If direct reaction check fails, look for pattern of OH → ether conversion
                if not is_ether_formation:
                    # Check for alcohol or phenol in reactants
                    has_oh_group = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Phenol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            has_oh_group = True
                            break

                    # Check for ether in product
                    has_ether = checker.check_fg("Ether", product)

                    # If both conditions are met, it's likely an ether formation
                    is_ether_formation = has_oh_group and has_ether

                if is_ether_formation:
                    ether_formations += 1
                    print(f"Ether formation detected at depth {depth}, reaction: {rsmi}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if multiple ether formations are detected
    return ether_formations >= 2
