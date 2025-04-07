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
    This function detects a 'build-couple-pair' strategy where a complex heterocycle is built
    before final functionalization steps.

    In retrosynthetic analysis:
    - Higher depth = earlier stage in forward synthesis (build phase)
    - Lower depth = later stage in forward synthesis (pair/functionalization phase)
    """
    # We need to track multiple aspects of the synthesis
    heterocycle_built = False
    coupling_phase = False
    late_functionalization = False
    heterocycle_build_depth = -1
    coupling_depth = -1
    functionalization_depth = -1

    # List of heterocycles to check for
    heterocycle_types = [
        "triazole",
        "tetrazole",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzotriazole",
    ]

    # List of functionalization reactions to check for
    functionalization_reactions = [
        "Sulfonamide synthesis (Schotten-Baumann) primary amine",
        "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "N-methylation",
        "O-methylation",
        "Schotten-Baumann to ester",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Esterification of Carboxylic Acids",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
    ]

    # List of coupling reactions to check for
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Sonogashira alkyne_aryl halide",
        "Heck terminal vinyl",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Stille reaction_aryl",
        "Negishi coupling",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_built, coupling_phase, late_functionalization
        nonlocal heterocycle_build_depth, coupling_depth, functionalization_depth

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for heterocycle formation (build phase)
                if not heterocycle_built:
                    product_has_heterocycle = False
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, product):
                            product_has_heterocycle = True
                            print(f"Found {heterocycle} in product at depth {depth}")
                            break

                    if product_has_heterocycle:
                        # Check if heterocycle was formed in this step (not present in reactants)
                        reactants_have_heterocycle = False
                        for reactant in reactants:
                            for heterocycle in heterocycle_types:
                                if checker.check_ring(heterocycle, reactant):
                                    reactants_have_heterocycle = True
                                    print(
                                        f"Found {heterocycle} in reactant at depth {depth}"
                                    )
                                    break
                            if reactants_have_heterocycle:
                                break

                        if not reactants_have_heterocycle:
                            heterocycle_built = True
                            heterocycle_build_depth = depth
                            print(f"Heterocycle built at depth {depth}")

                # Check for coupling phase (couple phase)
                # Check for coupling reactions at any point
                if not coupling_phase:
                    for reaction in coupling_reactions:
                        if checker.check_reaction(reaction, rsmi):
                            coupling_phase = True
                            coupling_depth = depth
                            print(
                                f"Coupling reaction detected at depth {depth}: {reaction}"
                            )
                            break

                # Check for late functionalization (pair phase)
                # In retrosynthetic analysis, lower depths are later stages in forward synthesis
                if not late_functionalization:
                    # Check for functionalization reactions
                    for reaction in functionalization_reactions:
                        if checker.check_reaction(reaction, rsmi):
                            late_functionalization = True
                            functionalization_depth = depth
                            print(
                                f"Functionalization detected at depth {depth}: {reaction}"
                            )
                            break

                    # Also check for specific functional groups being added
                    if not late_functionalization:
                        # Check for functional groups in product that aren't in reactants
                        fg_types = [
                            "Sulfonamide",
                            "Ester",
                            "Amide",
                            "Carbamate",
                            "Urea",
                            "Triflate",
                            "Mesylate",
                            "Tosylate",
                        ]
                        for fg in fg_types:
                            if checker.check_fg(fg, product):
                                # Check if this FG was added in this step
                                reactants_have_fg = False
                                for reactant in reactants:
                                    if checker.check_fg(fg, reactant):
                                        reactants_have_fg = True
                                        break

                                if not reactants_have_fg:
                                    late_functionalization = True
                                    functionalization_depth = depth
                                    print(
                                        f"Functionalization detected at depth {depth}: {fg} added"
                                    )
                                    break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # For a true build-couple-pair strategy:
    # 1. We need both heterocycle building and late functionalization
    # 2. In forward synthesis, heterocycle building should happen before functionalization
    # 3. In retrosynthetic analysis, heterocycle building depth should be greater than functionalization depth

    valid_strategy = False

    if heterocycle_built and late_functionalization:
        # Check if the heterocycle was built earlier in the forward synthesis (higher depth in retrosynthesis)
        # than the functionalization occurred
        if heterocycle_build_depth > functionalization_depth:
            valid_strategy = True
            print(
                f"Valid build-couple-pair strategy: heterocycle built at depth {heterocycle_build_depth}, functionalization at depth {functionalization_depth}"
            )
        else:
            print(
                f"Invalid sequence: heterocycle built at depth {heterocycle_build_depth}, functionalization at depth {functionalization_depth}"
            )

    # Coupling is optional but should be between build and pair phases if present
    if coupling_phase and valid_strategy:
        if functionalization_depth < coupling_depth < heterocycle_build_depth:
            print(f"Coupling phase detected at appropriate depth {coupling_depth}")
        else:
            print(
                f"Coupling detected at depth {coupling_depth}, but not between build and pair phases"
            )

    print(
        f"Final result: heterocycle_built={heterocycle_built}, coupling_phase={coupling_phase}, late_functionalization={late_functionalization}, valid_strategy={valid_strategy}"
    )
    return valid_strategy
