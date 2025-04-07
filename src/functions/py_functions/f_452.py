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
    Detects if the synthesis involves formation of an organotin reagent from an aryl halide.
    """
    found_organotin = False

    def dfs_traverse(node):
        nonlocal found_organotin

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction: {rsmi}")

                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if product contains tin
                if "Sn" in product_part:
                    print(f"Found product with tin: {product_part}")

                    # Check if any reactant has aromatic halide
                    for reactant in reactants:
                        print(f"Checking reactant: {reactant}")

                        # Check for aromatic halides
                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ):

                            # Verify this is a reaction that forms organotin compounds
                            if (
                                checker.check_reaction("Stille reaction_aryl", rsmi)
                                or checker.check_reaction("Stille reaction_vinyl", rsmi)
                                or checker.check_reaction(
                                    "Stille reaction_benzyl", rsmi
                                )
                                or checker.check_reaction("Stille reaction_allyl", rsmi)
                                or checker.check_reaction("Stille reaction_other", rsmi)
                            ):
                                print("Found organotin reagent formation from halide")
                                found_organotin = True
                                break

                            # If not a known Stille reaction, check if it's a general metal-halogen exchange
                            # that results in an organotin compound
                            product_mol = Chem.MolFromSmiles(product_part)
                            reactant_mol = Chem.MolFromSmiles(reactant)

                            if product_mol and reactant_mol:
                                # Check if the reaction involves replacing a halogen with tin
                                if (
                                    "Br" in reactant
                                    or "Cl" in reactant
                                    or "I" in reactant
                                ):
                                    print(
                                        "Found potential organotin reagent formation from halide"
                                    )
                                    found_organotin = True
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {found_organotin}")
    return found_organotin
