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
    This function detects if the synthesis involves indazole ring formation.
    """
    indazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal indazole_formation_detected

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if indazole is in product
                if checker.check_ring("indazole", product_smiles):
                    print(f"Indazole found in product: {product_smiles}")

                    # Check if indazole is not in any reactants
                    reactants_have_indazole = False
                    for reactant_smiles in reactants_smiles:
                        if checker.check_ring("indazole", reactant_smiles):
                            reactants_have_indazole = True
                            print(f"Indazole also found in reactant: {reactant_smiles}")
                            break

                    if not reactants_have_indazole:
                        # Check if this is a known indazole formation reaction
                        if checker.check_reaction(
                            "Intramolecular amination (heterocycle formation)", rsmi
                        ):
                            print(
                                f"Indazole formation via intramolecular amination detected: {rsmi}"
                            )
                            indazole_formation_detected = True
                        # Check for other potential indazole formation reactions
                        elif (
                            checker.check_reaction(
                                "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                            )
                            or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                            or checker.check_reaction(
                                "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                            )
                        ):
                            print(f"Indazole formation via cycloaddition detected: {rsmi}")
                            indazole_formation_detected = True
                        # If no specific reaction type is identified but indazole is formed
                        else:
                            print(
                                f"Indazole formation detected (unspecified reaction type): {rsmi}"
                            )
                            indazole_formation_detected = True
            except Exception as e:
                print(f"Error in indazole formation detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Indazole formation strategy result: {indazole_formation_detected}")
    return indazole_formation_detected
