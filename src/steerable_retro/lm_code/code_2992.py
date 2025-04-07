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
    Detects if the synthesis involves formation of a pyrazole ring from acyclic precursors.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if result:  # Optional early return if we already found a match
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if product contains pyrazole ring
                if checker.check_ring("pyrazole", product):
                    print(f"Found product with pyrazole ring: {product}")

                    # Check if any reactant contains hydrazine but not pyrazole
                    hydrazine_found = False
                    acyclic_precursor = False
                    carbonyl_found = False

                    for reactant in reactants:
                        # Check for hydrazine
                        if checker.check_fg("Hydrazine", reactant):
                            hydrazine_found = True
                            print(f"Found reactant with hydrazine: {reactant}")

                            # Ensure reactant doesn't already have pyrazole
                            if not checker.check_ring("pyrazole", reactant):
                                acyclic_precursor = True
                                print(f"Confirmed reactant is acyclic (no pyrazole): {reactant}")

                        # Check for carbonyl-containing compounds
                        if (
                            checker.check_fg("Aldehyde", reactant)
                            or checker.check_fg("Ketone", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Carboxylic acid", reactant)
                        ):
                            carbonyl_found = True
                            print(f"Found reactant with carbonyl group: {reactant}")

                    # Check if this is a pyrazole formation reaction
                    if hydrazine_found and acyclic_precursor:
                        # First check if it's a known pyrazole formation reaction
                        if checker.check_reaction("{pyrazole}", rsmi):
                            result = True
                            print(
                                "Detected pyrazole formation reaction using named reaction pattern"
                            )
                        # Alternative check for pyrazole reaction
                        elif carbonyl_found:
                            result = True
                            print(
                                "Detected pyrazole formation from hydrazine and carbonyl compound"
                            )
                        else:
                            print(
                                "Hydrazine and acyclic precursor found, but couldn't confirm pyrazole formation mechanism"
                            )
            except Exception as e:
                print(f"Error processing SMILES in pyrazole_formation_strategy: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
