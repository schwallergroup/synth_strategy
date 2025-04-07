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
    This function specifically detects the formation of an imidazopyridine scaffold.
    Imidazopyridine is a bicyclic heterocycle consisting of imidazole fused to pyridine.
    """
    imidazopyridine_formation_detected = False

    def dfs_traverse(node):
        nonlocal imidazopyridine_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for imidazopyridine in product
            # Imidazopyridine has both imidazole and pyridine rings fused together
            product_has_imidazole = checker.check_ring("imidazole", product_smiles)
            product_has_pyridine = checker.check_ring("pyridine", product_smiles)

            if product_has_imidazole and product_has_pyridine:
                print(f"Product contains imidazole and pyridine rings: {product_smiles}")

                # Check if all reactants have both rings
                all_reactants_have_both_rings = True
                for reactant in reactants_smiles:
                    if not (
                        checker.check_ring("imidazole", reactant)
                        and checker.check_ring("pyridine", reactant)
                    ):
                        all_reactants_have_both_rings = False
                        break

                # If not all reactants have both rings, this reaction likely formed the scaffold
                if not all_reactants_have_both_rings:
                    print(f"Potential imidazopyridine formation detected in reaction: {rsmi}")

                    # Check for specific reactions known to form imidazopyridines
                    if (
                        checker.check_reaction("benzimidazole_derivatives_aldehyde", rsmi)
                        or checker.check_reaction(
                            "benzimidazole_derivatives_carboxylic-acid/ester", rsmi
                        )
                        or checker.check_reaction("Pictet-Spengler", rsmi)
                        or checker.check_reaction("imidazole", rsmi)
                    ):
                        print(f"Confirmed imidazopyridine-forming reaction type")
                        imidazopyridine_formation_detected = True
                    else:
                        # Additional check: verify that at least one reactant has pyridine
                        # This helps confirm we're fusing an imidazole to an existing pyridine
                        reactant_has_pyridine = any(
                            checker.check_ring("pyridine", r) for r in reactants_smiles
                        )
                        if reactant_has_pyridine:
                            print(
                                f"Generic imidazopyridine formation detected (pyridine + imidazole fusion)"
                            )
                            imidazopyridine_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return imidazopyridine_formation_detected
