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
    Detects a synthetic strategy involving late-stage N-alkylation of a piperazine or similar nitrogen heterocycle.
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            is_late_stage = depth <= 1

            if is_late_stage:
                print(f"Checking late-stage reaction at depth {depth}: {rsmi}")

                # Check if this is an N-alkylation reaction
                is_n_alkylation = (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction(
                        "Eschweiler-Clarke Primary Amine Methylation", rsmi
                    )
                    or checker.check_reaction(
                        "Eschweiler-Clarke Secondary Amine Methylation", rsmi
                    )
                    or checker.check_reaction("N-methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                )  # Sometimes used for N-alkylation too

                print(f"Is N-alkylation reaction: {is_n_alkylation}")

                # If not detected by reaction checkers, try to detect by examining reactants and products
                if not is_n_alkylation:
                    # Check if product has a tertiary amine that wasn't in reactants
                    has_tertiary_amine_product = checker.check_fg(
                        "Tertiary amine", product
                    )

                    if has_tertiary_amine_product:
                        # Check if any reactant has a secondary amine
                        has_secondary_amine_reactant = any(
                            checker.check_fg("Secondary amine", r) for r in reactants
                        )

                        if has_secondary_amine_reactant:
                            print(
                                f"Detected N-alkylation pattern: secondary amine → tertiary amine"
                            )
                            is_n_alkylation = True

                if is_n_alkylation:
                    print(f"N-alkylation reaction detected: {rsmi}")

                    # Check if the product contains a nitrogen heterocycle
                    nitrogen_heterocycles = [
                        "piperazine",
                        "piperidine",
                        "pyrrolidine",
                        "morpholine",
                        "azepane",
                        "azetidine",
                        "pyridine",
                        "pyrrole",
                        "imidazole",
                        "triazole",
                        "tetrazole",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "purine",
                        "carbazole",
                        "acridine",
                        "benzimidazole",
                        "benzotriazole",
                        "indazole",
                    ]

                    # Check product for nitrogen heterocycles
                    product_heterocycles = []
                    for ring in nitrogen_heterocycles:
                        if checker.check_ring(ring, product):
                            product_heterocycles.append(ring)

                    if product_heterocycles:
                        print(
                            f"Nitrogen heterocycles found in product: {product_heterocycles}"
                        )

                        # Check reactants for the same heterocycles
                        for reactant in reactants:
                            reactant_heterocycles = []
                            for ring in product_heterocycles:
                                if checker.check_ring(ring, reactant):
                                    reactant_heterocycles.append(ring)

                            if reactant_heterocycles:
                                print(
                                    f"Nitrogen heterocycles found in reactant: {reactant_heterocycles}"
                                )

                                # Check if the reactant has a secondary amine in the heterocycle
                                if checker.check_fg("Secondary amine", reactant):
                                    print(
                                        f"Secondary amine found in reactant with heterocycle"
                                    )
                                    n_alkylation_detected = True
                                    break
                                # Check if the reactant has a primary amine in the heterocycle
                                elif checker.check_fg("Primary amine", reactant):
                                    print(
                                        f"Primary amine found in reactant with heterocycle"
                                    )
                                    n_alkylation_detected = True
                                    break

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"N-alkylation strategy detected: {n_alkylation_detected}")
    return n_alkylation_detected
