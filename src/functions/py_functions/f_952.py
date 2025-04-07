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
    Detects if the synthesis route includes a late-stage ester hydrolysis
    (conversion of ester to carboxylic acid in the final step)
    """
    ester_hydrolysis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_found

        # Set depth for the current node
        node["depth"] = depth

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for various ester hydrolysis reaction types
                is_hydrolysis = checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    rsmi,
                )
                is_saponification_methyl = checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                )
                is_saponification_alkyl = checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                )

                if is_hydrolysis or is_saponification_methyl or is_saponification_alkyl:
                    reaction_type = "hydrolysis" if is_hydrolysis else "saponification"
                    print(f"Found {reaction_type} reaction: {rsmi}")

                    # Verify ester in reactants and carboxylic acid in product
                    ester_in_reactants = False
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            ester_in_reactants = True
                            print(f"Found ester in reactant: {reactant}")
                            break

                    acid_in_product = checker.check_fg("Carboxylic acid", product)
                    if acid_in_product:
                        print(f"Found carboxylic acid in product: {product}")

                    if ester_in_reactants and acid_in_product:
                        print("Confirmed late-stage ester hydrolysis")
                        ester_hydrolysis_found = True
                else:
                    # Manual check for ester to acid conversion if reaction type check failed
                    print(
                        "Standard reaction checks failed, performing manual verification"
                    )
                    ester_in_reactants = any(
                        checker.check_fg("Ester", reactant) for reactant in reactants
                    )
                    acid_in_product = checker.check_fg("Carboxylic acid", product)

                    # Additional check: look for sodium or other bases in reagents
                    reagents = rsmi.split(">")[1].split(".")
                    has_base = any("[Na+]" in reagent for reagent in reagents)
                    has_water = any(reagent.strip() == "O" for reagent in reagents)

                    if (
                        ester_in_reactants
                        and acid_in_product
                        and (has_base or has_water)
                    ):
                        print(
                            "Manually confirmed ester hydrolysis with base/water present"
                        )
                        ester_hydrolysis_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return ester_hydrolysis_found
