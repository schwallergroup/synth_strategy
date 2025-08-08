#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy involving late-stage amide coupling
    with a heterocyclic partner.
    """
    print("Starting late_stage_amide_coupling_strategy analysis")
    late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_coupling

        # Check only late-stage reactions (depth 0, 1, or 2)
        if (
            node["type"] == "reaction"
            and depth <= 2
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide coupling reaction using the checker
            is_amide_coupling = any(
                [
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    ),
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    ),
                    checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi),
                    checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    ),
                    checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi),
                    checker.check_reaction("Ester with primary amine to amide", rsmi),
                    checker.check_reaction("Ester with secondary amine to amide", rsmi),
                    checker.check_reaction("Schotten-Baumann_amide", rsmi),
                    checker.check_reaction("Acylation of primary amines", rsmi),
                    checker.check_reaction("Acylation of secondary amines", rsmi),
                ]
            )

            # If no specific reaction type is detected, check for general amide formation pattern
            if not is_amide_coupling:
                # Check if product has amide that wasn't in reactants
                has_amide_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                has_amide_reactants = any(
                    [
                        (
                            checker.check_fg("Primary amide", r)
                            or checker.check_fg("Secondary amide", r)
                            or checker.check_fg("Tertiary amide", r)
                        )
                        for r in reactants
                    ]
                )

                # Check for acid and amine components
                has_acid = any([checker.check_fg("Carboxylic acid", r) for r in reactants])
                has_amine = any(
                    [
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    ]
                )

                # If product has amide, reactants have acid and amine, and reactants don't have amide
                if has_amide_product and has_acid and has_amine and not has_amide_reactants:
                    print("Detected amide formation based on functional group analysis")
                    is_amide_coupling = True

            if is_amide_coupling:
                print(f"Found amide coupling reaction at depth {depth}: {rsmi}")

                # Check for amide in product
                has_amide_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                if has_amide_product:
                    print(f"Product contains amide functional group: {product}")

                    # Check for carboxylic acid or activated derivatives in reactants
                    acid_reactants = []
                    for i, r in enumerate(reactants):
                        if (
                            checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Acyl halide", r)
                            or checker.check_fg("Ester", r)
                        ):
                            acid_reactants.append((i, r))
                            print(f"Found acid component in reactant {i}: {r}")

                    # Check for amine in reactants
                    amine_reactants = []
                    for i, r in enumerate(reactants):
                        if (
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                        ):
                            amine_reactants.append((i, r))
                            print(f"Found amine component in reactant {i}: {r}")

                    # Define heterocycle types to check
                    heterocycle_types = [
                        "pyridine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrazole",
                        "triazole",
                        "tetrazole",
                        "pyrimidine",
                        "pyrazine",
                        "indole",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "quinoline",
                        "isoquinoline",
                        "piperidine",
                        "piperazine",
                        "morpholine",
                        "thiomorpholine",
                    ]

                    # Check for heterocycles in reactants and track which ones
                    heterocycle_in_reactants = {}
                    for i, r in enumerate(reactants):
                        for ring_type in heterocycle_types:
                            if checker.check_ring(ring_type, r):
                                if i not in heterocycle_in_reactants:
                                    heterocycle_in_reactants[i] = []
                                heterocycle_in_reactants[i].append(ring_type)
                                print(f"Found heterocycle ({ring_type}) in reactant {i}: {r}")

                    # Check if heterocycle is in the product as well
                    heterocycle_in_product = []
                    for ring_type in heterocycle_types:
                        if checker.check_ring(ring_type, product):
                            heterocycle_in_product.append(ring_type)
                            print(f"Found heterocycle ({ring_type}) in product: {product}")

                    # Determine if we have a valid late-stage amide coupling with heterocyclic partner
                    has_acid_derivative = len(acid_reactants) > 0
                    has_amine = len(amine_reactants) > 0

                    print(f"Heterocycle reactants: {heterocycle_in_reactants}")
                    print(f"Acid reactants: {acid_reactants}")
                    print(f"Amine reactants: {amine_reactants}")

                    # Check if any heterocycle is in either the acid or amine component
                    heterocycle_in_coupling = False
                    for reactant_idx in heterocycle_in_reactants:
                        for acid_idx, _ in acid_reactants:
                            if reactant_idx == acid_idx:
                                heterocycle_in_coupling = True
                                print(
                                    f"Heterocycle found in acid component (reactant {reactant_idx})"
                                )
                        for amine_idx, _ in amine_reactants:
                            if reactant_idx == amine_idx:
                                heterocycle_in_coupling = True
                                print(
                                    f"Heterocycle found in amine component (reactant {reactant_idx})"
                                )

                    # Also consider heterocycle in product if it wasn't in reactants
                    if not heterocycle_in_coupling and heterocycle_in_product:
                        heterocycle_in_coupling = True
                        print("Heterocycle found in product but not in specific reactants")

                    if has_acid_derivative and has_amine and heterocycle_in_coupling:
                        print(
                            f"✓ Confirmed late-stage amide coupling with heterocyclic partner at depth {depth}"
                        )
                        late_stage_amide_coupling = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Final result: late_stage_amide_coupling = {late_stage_amide_coupling}")
    return late_stage_amide_coupling
