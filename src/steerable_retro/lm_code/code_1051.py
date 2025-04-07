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
    This function detects if the synthetic route involves oxazole ring formation.
    """
    oxazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal oxazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if oxazole is in product
                if checker.check_ring("oxazole", product):
                    print(f"Found oxazole in product: {product}")

                    # Check if oxazole is not in any reactant
                    reactant_has_oxazole = False
                    for reactant in reactants:
                        if checker.check_ring("oxazole", reactant):
                            print(f"Found oxazole in reactant: {reactant}")
                            reactant_has_oxazole = True
                            break

                    if not reactant_has_oxazole:
                        # Check if this is a known oxazole formation reaction
                        if (
                            checker.check_reaction("benzoxazole formation from aldehyde", rsmi)
                            or checker.check_reaction(
                                "benzoxazole formation from acyl halide", rsmi
                            )
                            or checker.check_reaction(
                                "benzoxazole formation from ester/carboxylic acid", rsmi
                            )
                            or checker.check_reaction(
                                "benzoxazole formation (intramolecular)", rsmi
                            )
                        ):
                            print(f"Oxazole formation reaction detected: {rsmi}")
                            oxazole_formation_detected = True
                        else:
                            # Generic check for oxazole formation
                            print(f"Potential oxazole formation detected: {rsmi}")
                            oxazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Oxazole formation detected: {oxazole_formation_detected}")
    return oxazole_formation_detected
