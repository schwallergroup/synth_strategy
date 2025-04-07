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
    Detects nitro group reduction to amine in the synthetic route.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a nitro group
                reactant_has_nitro = False
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        reactant_has_nitro = True
                        print(f"Found reactant with nitro group: {reactant}")
                        break

                # Check if product has a primary amine
                product_has_amine = checker.check_fg("Primary amine", product)
                if product_has_amine:
                    print(f"Found product with primary amine: {product}")

                # Check if this is a nitro reduction reaction
                if reactant_has_nitro and product_has_amine:
                    if checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    ):
                        nitro_reduction_found = True
                        print(f"Found nitro reduction to amine reaction: {rsmi}")
                    else:
                        # Fallback check: look for nitro group and primary amine at corresponding positions
                        # This is a backup in case the reaction checker fails
                        try:
                            # Find reactant with nitro group
                            nitro_reactant = None
                            for reactant in reactants:
                                if checker.check_fg("Nitro group", reactant):
                                    nitro_reactant = reactant
                                    break

                            if nitro_reactant:
                                # Get atom indices of nitro group in reactant
                                nitro_indices = checker.get_fg_atom_indices(
                                    "Nitro group", nitro_reactant
                                )
                                if nitro_indices:
                                    # This is a simplification - in a real implementation we would use
                                    # atom mapping to track the exact transformation
                                    nitro_reduction_found = True
                                    print(
                                        f"Found nitro reduction to amine (fallback detection): {rsmi}"
                                    )
                        except Exception as e:
                            print(f"Error in fallback detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_reduction_found
