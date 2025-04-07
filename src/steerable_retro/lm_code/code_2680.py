#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects Sonogashira coupling (aryl/vinyl halide/triflate + terminal alkyne) in the synthesis.
    """
    sonogashira_found = False

    def dfs_traverse(node):
        nonlocal sonogashira_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for Sonogashira coupling using the checker function
            sonogashira_reaction_types = [
                "Sonogashira acetylene_aryl halide",
                "Sonogashira alkyne_aryl halide",
                "Sonogashira acetylene_aryl OTf",
                "Sonogashira alkyne_aryl OTf",
                "Sonogashira acetylene_alkenyl halide",
                "Sonogashira alkyne_alkenyl halide",
                "Sonogashira acetylene_alkenyl OTf",
                "Sonogashira alkyne_alkenyl OTf",
                "Sonogashira acetylene_acyl halide",
                "Sonogashira alkyne_acyl halide",
            ]

            for reaction_type in sonogashira_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Sonogashira coupling detected: {reaction_type}")
                    sonogashira_found = True
                    break

            # If no direct reaction match, check for the characteristic pattern
            if not sonogashira_found:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if len(reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for aryl-alkyne bond in product
                    if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c-C#C")):
                        # Check for aryl/vinyl halide/triflate and terminal alkyne in reactants
                        electrophile_found = False
                        alkyne_found = False

                        for reactant in reactants:
                            if (
                                checker.check_fg("Aromatic halide", reactant)
                                or checker.check_fg("Alkenyl halide", reactant)
                                or checker.check_fg("Triflate", reactant)
                            ):
                                electrophile_found = True
                                print(f"Electrophile found in reactant: {reactant}")

                            if checker.check_fg("Alkyne", reactant):
                                alkyne_found = True
                                print(f"Alkyne found in reactant: {reactant}")

                        if electrophile_found and alkyne_found:
                            sonogashira_found = True
                            print("Sonogashira coupling pattern detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return sonogashira_found
