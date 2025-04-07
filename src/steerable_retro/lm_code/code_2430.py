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
    This function detects if the synthesis route involves nucleophilic aromatic substitution,
    specifically C-Cl to C-N transformation on a pyridine ring.
    """
    nas_detected = False

    def dfs_traverse(node):
        nonlocal nas_detected

        if node["type"] == "reaction" and not nas_detected:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # First check if this is a nucleophilic substitution reaction
                if (
                    checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                ):
                    print(f"Detected nucleophilic substitution reaction: {rsmi}")

                    # Check for pyridine ring in reactants
                    pyridine_in_reactants = any(
                        checker.check_ring("pyridine", r) for r in reactants_smiles
                    )

                    # Check for aromatic halide in reactants
                    aromatic_halide_in_reactants = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    )

                    # Check for amine-containing product
                    amine_in_product = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    if pyridine_in_reactants and aromatic_halide_in_reactants and amine_in_product:
                        print(
                            "Confirmed: Nucleophilic aromatic substitution on pyridine (C-Cl to C-N)"
                        )
                        nas_detected = True

                # If not detected by reaction type, try to infer from functional group changes
                if not nas_detected:
                    # Check for pyridine ring and aromatic halide in reactants
                    for reactant_smiles in reactants_smiles:
                        if checker.check_ring("pyridine", reactant_smiles) and checker.check_fg(
                            "Aromatic halide", reactant_smiles
                        ):
                            # Check if product has pyridine ring and amine group
                            if checker.check_ring("pyridine", product_smiles) and (
                                checker.check_fg("Primary amine", product_smiles)
                                or checker.check_fg("Secondary amine", product_smiles)
                                or checker.check_fg("Tertiary amine", product_smiles)
                                or checker.check_fg("Aniline", product_smiles)
                            ):
                                print(
                                    f"Inferred nucleophilic aromatic substitution from functional group changes: {rsmi}"
                                )
                                nas_detected = True
                                break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    return nas_detected
