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
    This function detects a late-stage nucleophilic aromatic substitution with an amine.
    """
    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar

        print(f"Examining node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Late stage (expanded depth check)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction SMILES: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr or equivalent reactions with amines
                is_snar = (
                    checker.check_reaction("Nucleophilic substitution aromatic", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution amine", rsmi
                    )
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction(
                        "Goldberg coupling aryl amine-aryl chloride", rsmi
                    )
                    or checker.check_reaction(
                        "Goldberg coupling aryl amide-aryl chloride", rsmi
                    )
                    or checker.check_reaction("Ullmann condensation", rsmi)
                )

                print(f"Is SNAr or related reaction: {is_snar}")

                if is_snar:
                    print(
                        f"Found potential SNAr or N-arylation reaction at depth {depth}"
                    )

                    # Track reactants with specific functional groups
                    aromatic_halide_reactant = None
                    amine_reactant = None

                    for reactant in reactants:
                        print(f"Checking reactant: {reactant}")

                        if checker.check_fg("Aromatic halide", reactant):
                            aromatic_halide_reactant = reactant
                            print(f"Found aromatic halide in reactant: {reactant}")

                        # Check for any type of amine or nitrogen nucleophile
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            amine_reactant = reactant
                            print(f"Found amine in reactant: {reactant}")

                    # Verify both required reactants are present
                    if aromatic_halide_reactant and amine_reactant:
                        # Verify the product no longer has the aromatic halide
                        if not checker.check_fg("Aromatic halide", product):
                            # The reaction has been identified as SNAr with amine
                            print(f"Confirmed SNAr with amine at depth {depth}")
                            found_snar = True
                else:
                    # Try a more general approach if specific reaction checks fail
                    print("Trying general approach for SNAr detection")
                    aromatic_halide_found = False
                    amine_found = False

                    for reactant in reactants:
                        if checker.check_fg("Aromatic halide", reactant):
                            aromatic_halide_found = True
                            print(f"Found aromatic halide in reactant: {reactant}")

                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            amine_found = True
                            print(f"Found amine in reactant: {reactant}")

                    if aromatic_halide_found and amine_found:
                        # Check if product has lost the aromatic halide
                        if not checker.check_fg("Aromatic halide", product):
                            print(
                                f"Detected potential SNAr reaction by functional group analysis at depth {depth}"
                            )
                            found_snar = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_snar
