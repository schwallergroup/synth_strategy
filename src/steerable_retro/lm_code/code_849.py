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
    This function detects alcohol protection as ether in the synthetic route.
    """
    alcohol_protection_detected = False

    def dfs_traverse(node):
        nonlocal alcohol_protection_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if this is a known protection reaction
                if checker.check_reaction(
                    "Williamson Ether Synthesis", rsmi
                ) or checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    print("Detected alcohol protection via known reaction")
                    alcohol_protection_detected = True
                    return

                # Check for alcohol in reactants
                alcohol_in_reactants = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                )

                # Check for ether in product and absence of alcohol
                ether_in_product = checker.check_fg("Ether", product)
                alcohol_in_product = (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Aromatic alcohol", product)
                )

                # Check for halide or other leaving groups in reactants
                leaving_group_present = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    or checker.check_fg("Triflate", r)
                    or checker.check_fg("Mesylate", r)
                    or checker.check_fg("Tosylate", r)
                    for r in reactants
                )

                if (
                    alcohol_in_reactants
                    and ether_in_product
                    and (not alcohol_in_product or alcohol_in_reactants != alcohol_in_product)
                    and leaving_group_present
                ):
                    print("Detected alcohol protection as ether via functional group analysis")
                    alcohol_protection_detected = True
            except Exception as e:
                print(f"Error in alcohol protection detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return alcohol_protection_detected
