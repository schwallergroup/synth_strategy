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
    This function detects an amide coupling strategy in the synthesis.
    It looks for the formation of an amide bond (C(=O)N) from an amine and
    a carboxylic acid or derivative.
    """
    amide_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction directly
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    )
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Ester with secondary amine to amide", rsmi
                    )
                ):
                    print(f"Found amide coupling reaction at depth {depth}: {rsmi}")
                    amide_coupling_detected = True
                else:
                    # Check for amine and acid/derivative in reactants and amide in product
                    amine_found = False
                    acid_derivative_found = False

                    for reactant in reactants:
                        if checker.check_fg(
                            "Primary amine", reactant
                        ) or checker.check_fg("Secondary amine", reactant):
                            amine_found = True
                            print(f"Found amine in reactant: {reactant}")

                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            acid_derivative_found = True
                            print(f"Found acid derivative in reactant: {reactant}")

                    # Check for amide in product
                    if amine_found and acid_derivative_found:
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            print(f"Found amide in product: {product}")

                            # Verify amide wasn't already present in all reactants
                            all_reactants_have_amide = True
                            for reactant in reactants:
                                if not (
                                    checker.check_fg("Primary amide", reactant)
                                    or checker.check_fg("Secondary amide", reactant)
                                    or checker.check_fg("Tertiary amide", reactant)
                                ):
                                    all_reactants_have_amide = False
                                    break

                            if not all_reactants_have_amide:
                                print(f"Found amide coupling at depth {depth}")
                                amide_coupling_detected = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return amide_coupling_detected
