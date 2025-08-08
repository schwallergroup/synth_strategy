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
    This function detects if the synthesis involves a Sonogashira-type disconnection
    (aryl/vinyl halide + terminal alkyne).
    """
    has_sonogashira_disconnection = False

    def dfs_traverse(node):
        nonlocal has_sonogashira_disconnection

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for any Sonogashira reaction type
                if (
                    checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_aryl OTf", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_alkenyl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_alkenyl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_alkenyl OTf", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_alkenyl OTf", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_acyl halide", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_acyl halide", rsmi)
                ):

                    print(f"Found Sonogashira reaction: {rsmi}")
                    has_sonogashira_disconnection = True

                # If no direct Sonogashira reaction is detected, check for the characteristic functional groups
                if not has_sonogashira_disconnection:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    if len(reactants_smiles) > 1:
                        # Check if one reactant has aryl/vinyl halide and another has terminal alkyne
                        has_aryl_halide = any(
                            checker.check_fg("Aromatic halide", r)
                            or checker.check_fg("Alkenyl halide", r)
                            for r in reactants_smiles
                        )
                        has_terminal_alkyne = any(
                            checker.check_fg("Alkyne", r) for r in reactants_smiles
                        )
                        has_internal_alkyne = checker.check_fg("Alkyne", product_smiles)

                        if has_aryl_halide and has_terminal_alkyne and has_internal_alkyne:
                            print(f"Found Sonogashira-like disconnection: {rsmi}")
                            has_sonogashira_disconnection = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Sonogashira disconnection strategy detected: {has_sonogashira_disconnection}")
    return has_sonogashira_disconnection
