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
    This function detects Suzuki coupling (aryl halide + boronic acid/ester) in the synthetic route.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if suzuki_coupling_found:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a Suzuki coupling reaction using the checker
            suzuki_reaction_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki",
            ]

            # First check if the reaction is explicitly a Suzuki coupling
            for reaction_type in suzuki_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Suzuki coupling detected: {reaction_type} - {rsmi}")
                    suzuki_coupling_found = True
                    return

            # If not explicitly labeled, check if it has the characteristic functional groups
            try:
                reactants = rsmi.split(">")[0].split(".")

                # Check for required functional groups in reactants
                has_aryl_halide = False
                has_boronic = False

                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant):
                        has_aryl_halide = True
                        print(f"Found aromatic halide in reactant: {reactant}")

                    if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                        "Boronic ester", reactant
                    ):
                        has_boronic = True
                        print(f"Found boronic acid/ester in reactant: {reactant}")

                # If both required functional groups are present, it's likely a Suzuki coupling
                if has_aryl_halide and has_boronic:
                    # Double-check if it's a C-C bond formation reaction
                    print(
                        f"Detected unlabeled Suzuki coupling with required functional groups: {rsmi}"
                    )
                    suzuki_coupling_found = True
                    return
            except Exception as e:
                print(f"Error analyzing reactants: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Suzuki coupling found: {suzuki_coupling_found}")
    return suzuki_coupling_found
