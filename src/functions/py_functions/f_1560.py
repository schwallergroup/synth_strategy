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
    This function detects if the synthesis route involves a late-stage esterification
    (carboxylic acid to ester transformation in the final or penultimate step).

    In retrosynthesis, we're looking for ester → carboxylic acid transformations
    since we're going backwards from the target molecule.
    """
    esterification_found = False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, we check for ester in reactants and acid in product
                reactant_has_ester = any(
                    checker.check_fg("Ester", r) for r in reactants
                )
                product_has_acid = checker.check_fg("Carboxylic acid", product)

                print(f"  Reactant has ester: {reactant_has_ester}")
                print(f"  Product has carboxylic acid: {product_has_acid}")

                # Check if this is an esterification reaction (in forward direction)
                is_esterification = (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidative esterification of primary alcohols", rsmi
                    )
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction(
                        "Acetic anhydride and alcohol to ester", rsmi
                    )
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                )

                # Check if this is a hydrolysis reaction (esterification in reverse)
                is_hydrolysis = (
                    checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Ester saponification (alkyl deprotection)", rsmi
                    )
                )

                print(f"  Is esterification reaction: {is_esterification}")
                print(f"  Is hydrolysis reaction: {is_hydrolysis}")

                # In retrosynthesis, we're looking for ester → acid transformation
                # This can be either a forward esterification or a reverse hydrolysis
                if (
                    reactant_has_ester
                    and product_has_acid
                    and (is_esterification or is_hydrolysis)
                ):
                    print(f"Late-stage esterification detected at depth {depth}")
                    esterification_found = True

                # If we have the right functional groups but reaction type isn't detected,
                # do a more general check for ester to acid transformation
                elif reactant_has_ester and product_has_acid:
                    print(
                        f"Potential ester to acid transformation detected at depth {depth}"
                    )
                    # This is a fallback in case the specific reaction types aren't detected
                    esterification_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal to find late-stage esterification...")
    dfs_traverse(route)
    print(f"Esterification found: {esterification_found}")
    return esterification_found
