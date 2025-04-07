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
    This function detects a convergent synthesis where two independently
    functionalized fragments are combined in a late-stage coupling.
    """
    # Track if we found the pattern
    found_convergent_pattern = False
    # Track if we found the required transformations
    found_oxidation = False
    found_chlorination = False
    found_cross_coupling = False
    found_nucleophilic_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_pattern, found_oxidation, found_chlorination
        nonlocal found_cross_coupling, found_nucleophilic_substitution

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is the final coupling reaction (depth 1)
            if depth == 1:
                print("Checking for late-stage coupling...")
                # Check if we have two complex fragments being combined
                if len(reactants) >= 2:
                    # Check for methoxyphenyl and pyrazolopyrimidine fragments
                    has_methoxy = False
                    has_pyrazolo = False

                    for reactant in reactants:
                        # Check for methoxyphenyl fragment
                        if checker.check_fg("Ether", reactant) and checker.check_ring(
                            "benzene", reactant
                        ):
                            print(f"Found methoxyphenyl fragment in {reactant}")
                            has_methoxy = True

                        # Check for pyrazolopyrimidine-like structure
                        if checker.check_ring("pyrazole", reactant) or checker.check_ring(
                            "pyrimidine", reactant
                        ):
                            print(f"Found pyrazolopyrimidine fragment in {reactant}")
                            has_pyrazolo = True

                    # Check if the product contains both fragments
                    if has_methoxy and has_pyrazolo:
                        print(
                            "Found late-stage coupling of methoxyphenyl and pyrazolopyrimidine fragments"
                        )
                        found_convergent_pattern = True

                    # Check for nucleophilic substitution at coupling step
                    if (
                        checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                        or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                        or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    ):
                        print("Found nucleophilic substitution reaction at coupling step")
                        found_nucleophilic_substitution = True
                    else:
                        # Broader check for nucleophilic substitution
                        for reactant in reactants:
                            if (
                                "Cl" in reactant
                                or "Br" in reactant
                                or "I" in reactant
                                or checker.check_fg("Aromatic halide", reactant)
                                or checker.check_fg("Primary halide", reactant)
                            ):
                                # If we have a halide in a reactant and the product combines fragments
                                if has_methoxy and has_pyrazolo:
                                    print(
                                        "Found potential nucleophilic substitution at coupling step"
                                    )
                                    found_nucleophilic_substitution = True

            # Check for oxidation reactions at any depth
            if depth >= 2 and depth <= 5:
                print(f"Checking for oxidation at depth {depth}...")
                if (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    )
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                ):
                    print("Found alcohol oxidation reaction")
                    found_oxidation = True
                else:
                    # Backup check using functional groups
                    alcohol_in_reactant = False
                    carbonyl_in_product = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            alcohol_in_reactant = True
                            break

                    if (
                        checker.check_fg("Ketone", product)
                        or checker.check_fg("Aldehyde", product)
                        or checker.check_fg("Carboxylic acid", product)
                    ):
                        carbonyl_in_product = True

                    if alcohol_in_reactant and carbonyl_in_product:
                        print(f"Found alcohol to carbonyl oxidation at depth {depth}")
                        found_oxidation = True

            # Check for chlorination at any depth
            if depth >= 1 and depth <= 5:
                print(f"Checking for chlorination at depth {depth}...")
                if checker.check_reaction("Chlorination", rsmi) or checker.check_reaction(
                    "Aromatic chlorination", rsmi
                ):
                    print("Found chlorination reaction")
                    found_chlorination = True
                else:
                    # Check for chlorine in the product and reactants
                    has_chlorine_product = "Cl" in product

                    # Check for chlorinating reagents
                    chlorinating_reagent = False
                    for reactant in reactants:
                        if (
                            "N1C(=O)CCC(=O)C1Cl" in reactant
                            or "N-chloro" in reactant.lower()
                            or "ClCCl" in reactant
                        ):
                            chlorinating_reagent = True
                            break

                    if has_chlorine_product and chlorinating_reagent:
                        print(f"Found chlorination at depth {depth}")
                        found_chlorination = True

                # Check for chlorine in the pyrazolopyrimidine fragment
                for reactant in reactants:
                    if checker.check_ring("pyrazole", reactant) or checker.check_ring(
                        "pyrimidine", reactant
                    ):
                        if "Cl" in reactant:
                            print(f"Found chlorinated pyrazolopyrimidine at depth {depth}")
                            found_chlorination = True

            # Check for cross-coupling (Suzuki) at any depth
            if depth >= 2 and depth <= 5:
                print(f"Checking for Suzuki coupling at depth {depth}...")
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                ):
                    print("Found Suzuki cross-coupling reaction")
                    found_cross_coupling = True
                else:
                    # Backup check using functional groups
                    has_halide = False
                    has_boronic = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or "Br" in reactant
                            or "Cl" in reactant
                        ):
                            print(f"Found aromatic halide in {reactant}")
                            has_halide = True
                        if (
                            checker.check_fg("Boronic acid", reactant)
                            or checker.check_fg("Boronic ester", reactant)
                            or "B(O)" in reactant
                        ):
                            print(f"Found boronic acid/ester in {reactant}")
                            has_boronic = True

                    if has_halide and has_boronic:
                        print(f"Found cross-coupling reaction (Suzuki) at depth {depth}")
                        found_cross_coupling = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Pattern found: {found_convergent_pattern}")
    print(
        f"Transformations found: Oxidation={found_oxidation}, Chlorination={found_chlorination}, "
        f"Cross-coupling={found_cross_coupling}, Nucleophilic substitution={found_nucleophilic_substitution}"
    )

    # Return True if we found the convergent pattern and at least 3 of the 4 transformations
    transformations_found = sum(
        [found_oxidation, found_chlorination, found_cross_coupling, found_nucleophilic_substitution]
    )

    return found_convergent_pattern and transformations_found >= 3
