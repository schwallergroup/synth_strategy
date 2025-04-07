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

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis route involves a Suzuki coupling reaction.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check if the reaction is a Suzuki coupling using the checker function
            suzuki_reaction_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki",
            ]

            for reaction_type in suzuki_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Suzuki coupling detected at depth {depth}: {reaction_type}")
                    suzuki_detected = True
                    return  # No need to check further for this node

            # If direct reaction check didn't work, try checking for characteristic functional groups
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            has_boronic_acid = False
            has_boronic_ester = False
            has_aryl_halide = False

            for reactant in reactants:
                if checker.check_fg("Boronic acid", reactant):
                    has_boronic_acid = True
                    print(f"Found boronic acid in reactant: {reactant}")
                if checker.check_fg("Boronic ester", reactant):
                    has_boronic_ester = True
                    print(f"Found boronic ester in reactant: {reactant}")
                if checker.check_fg("Aromatic halide", reactant):
                    has_aryl_halide = True
                    print(f"Found aromatic halide in reactant: {reactant}")

            # If we have the characteristic reactants of a Suzuki coupling
            if (has_boronic_acid or has_boronic_ester) and has_aryl_halide:
                # Check if the product doesn't have the boronic acid/ester and halide anymore
                product_has_boronic_acid = checker.check_fg("Boronic acid", product)
                product_has_boronic_ester = checker.check_fg("Boronic ester", product)
                product_has_aryl_halide = checker.check_fg("Aromatic halide", product)

                if (
                    not (product_has_boronic_acid or product_has_boronic_ester)
                    and not product_has_aryl_halide
                ):
                    print(
                        f"Potential Suzuki coupling detected at depth {depth} based on reactants and product"
                    )
                    suzuki_detected = True
                    return  # No need to check further for this node

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)
            if suzuki_detected:
                return  # Early termination if we found what we're looking for

    dfs_traverse(route)
    return suzuki_detected
