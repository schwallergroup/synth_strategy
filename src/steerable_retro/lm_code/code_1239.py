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
    Detects if the synthesis route involves early heterocycle formation,
    specifically pyrazole formation from 1,3-dicarbonyl and hydrazine.
    Early stage means at high depth values (e.g., depth 4 or higher).
    """
    early_pyrazole_formation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal early_pyrazole_formation

        # Add depth to node metadata for debugging
        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = current_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is at a high depth (early in synthesis)
            if current_depth >= 4:
                # Check if this is a pyrazole formation reaction
                if checker.check_reaction("pyrazole", rsmi):
                    print(f"Detected pyrazole formation reaction at depth {current_depth}")
                    early_pyrazole_formation = True
                elif checker.check_ring("pyrazole", product):
                    print(f"Product contains pyrazole ring at depth {current_depth}")

                    # Check reactants for 1,3-dicarbonyl and hydrazine patterns
                    has_dicarbonyl = False
                    has_hydrazine = False

                    for reactant in reactants:
                        # Check for hydrazine or derivatives
                        if checker.check_fg("Hydrazine", reactant) or checker.check_fg(
                            "Hydrazone", reactant
                        ):
                            print(f"Found hydrazine/hydrazone in reactant: {reactant}")
                            has_hydrazine = True

                        # Check for 1,3-dicarbonyl pattern
                        # Since 1,3-dicarbonyl isn't in the FG list, we check for ketones/aldehydes
                        # and then verify the 1,3-dicarbonyl pattern
                        if checker.check_fg("Ketone", reactant) or checker.check_fg(
                            "Aldehyde", reactant
                        ):
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[#6](=[O])[#6][#6](=[O])")
                            ):
                                print(f"Found 1,3-dicarbonyl in reactant: {reactant}")
                                has_dicarbonyl = True

                    if has_dicarbonyl and has_hydrazine:
                        print(f"Confirmed early pyrazole formation at depth {current_depth}")
                        early_pyrazole_formation = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if early_pyrazole_formation:
        print("Early pyrazole formation strategy detected")
    else:
        print("No early pyrazole formation detected in the synthesis route")

    return early_pyrazole_formation
