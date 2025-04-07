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
    This function detects a sequential oxidation strategy where a vinyl group is
    transformed to an acid chloride via aldehyde and carboxylic acid intermediates,
    with the initial vinyl group introduced via Stille coupling.
    """
    # Track if we've found each transformation
    found_stille_coupling = False
    found_vinyl_to_aldehyde = False
    found_aldehyde_to_acid = False
    found_acid_to_acid_chloride = False

    # Track reaction depths for sequence verification
    stille_depth = -1
    vinyl_to_aldehyde_depth = -1
    aldehyde_to_acid_depth = -1
    acid_to_acid_chloride_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_stille_coupling, found_vinyl_to_aldehyde, found_aldehyde_to_acid, found_acid_to_acid_chloride
        nonlocal stille_depth, vinyl_to_aldehyde_depth, aldehyde_to_acid_depth, acid_to_acid_chloride_depth

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Stille coupling (introducing vinyl group)
            if checker.check_reaction("Stille reaction_vinyl", rsmi) or (
                any(checker.check_fg("Vinyl", r) for r in reactants)
                and any("Sn" in r for r in reactants)
                and checker.check_fg("Vinyl", product)
            ):
                print(f"Found Stille coupling at depth {depth}: {rsmi}")
                found_stille_coupling = True
                stille_depth = depth

            # Check for vinyl to aldehyde oxidation
            if (
                checker.check_reaction("Alkene oxidation to aldehyde", rsmi)
                or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
            ) or (
                any(checker.check_fg("Vinyl", r) for r in reactants)
                and checker.check_fg("Aldehyde", product)
                and not checker.check_fg("Carboxylic acid", product)
            ):
                print(f"Found vinyl to aldehyde oxidation at depth {depth}: {rsmi}")
                found_vinyl_to_aldehyde = True
                vinyl_to_aldehyde_depth = depth

            # Check for aldehyde to carboxylic acid oxidation
            if (
                checker.check_reaction(
                    "Oxidation of aldehydes to carboxylic acids", rsmi
                )
                or checker.check_reaction(
                    "Oxidation of aldehyde to carboxylic acid", rsmi
                )
            ) or (
                any(checker.check_fg("Aldehyde", r) for r in reactants)
                and checker.check_fg("Carboxylic acid", product)
            ):
                print(f"Found aldehyde to acid oxidation at depth {depth}: {rsmi}")
                found_aldehyde_to_acid = True
                aldehyde_to_acid_depth = depth

            # Check for acid to acid chloride conversion
            if any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            ) and checker.check_fg("Acyl halide", product):
                print(
                    f"Found acid to acid chloride conversion at depth {depth}: {rsmi}"
                )
                found_acid_to_acid_chloride = True
                acid_to_acid_chloride_depth = depth

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all transformations were found and in the correct sequence
    correct_sequence = (
        found_stille_coupling
        and found_vinyl_to_aldehyde
        and found_aldehyde_to_acid
        and found_acid_to_acid_chloride
        and stille_depth
        > vinyl_to_aldehyde_depth
        > aldehyde_to_acid_depth
        > acid_to_acid_chloride_depth
    )

    print(f"Stille coupling found: {found_stille_coupling} at depth {stille_depth}")
    print(
        f"Vinyl to aldehyde found: {found_vinyl_to_aldehyde} at depth {vinyl_to_aldehyde_depth}"
    )
    print(
        f"Aldehyde to acid found: {found_aldehyde_to_acid} at depth {aldehyde_to_acid_depth}"
    )
    print(
        f"Acid to acid chloride found: {found_acid_to_acid_chloride} at depth {acid_to_acid_chloride_depth}"
    )

    if correct_sequence:
        print("Found complete sequential vinyl to acid chloride strategy")
    else:
        print("Did not find complete sequential vinyl to acid chloride strategy")
        if not found_stille_coupling:
            print("Missing Stille coupling step")
        if not found_vinyl_to_aldehyde:
            print("Missing vinyl to aldehyde oxidation step")
        if not found_aldehyde_to_acid:
            print("Missing aldehyde to acid oxidation step")
        if not found_acid_to_acid_chloride:
            print("Missing acid to acid chloride conversion step")
        if not (
            stille_depth
            > vinyl_to_aldehyde_depth
            > aldehyde_to_acid_depth
            > acid_to_acid_chloride_depth
        ) and all(
            [
                found_stille_coupling,
                found_vinyl_to_aldehyde,
                found_aldehyde_to_acid,
                found_acid_to_acid_chloride,
            ]
        ):
            print("Steps found but not in correct sequence")

    return correct_sequence
