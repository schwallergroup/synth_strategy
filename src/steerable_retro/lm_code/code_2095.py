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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects the combined strategy of heterocycle-first synthesis
    with late-stage ether formation.
    """

    # Define helper functions for the individual strategies
    def heterocycle_first_strategy(route):
        """Check if heterocycles are formed in early stages of synthesis."""
        heterocycle_depths = []

        def dfs_heterocycle(node, depth=0):
            if node["type"] == "reaction":
                try:
                    rsmi = node["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]

                    # Check if the product contains any heterocyclic rings
                    heterocycle_rings = [
                        "furan",
                        "pyran",
                        "dioxane",
                        "tetrahydrofuran",
                        "tetrahydropyran",
                        "oxirane",
                        "oxetane",
                        "oxolane",
                        "oxane",
                        "dioxolane",
                        "dioxolene",
                        "trioxane",
                        "dioxepane",
                        "pyrrole",
                        "pyridine",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "pyrazine",
                        "pyridazine",
                        "triazole",
                        "tetrazole",
                        "pyrrolidine",
                        "piperidine",
                        "piperazine",
                        "morpholine",
                        "thiomorpholine",
                        "aziridine",
                        "azetidine",
                        "azepane",
                        "diazepane",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "purine",
                        "carbazole",
                        "acridine",
                        "thiophene",
                        "thiopyran",
                        "thiirane",
                        "thietane",
                        "thiolane",
                        "thiane",
                        "dithiane",
                        "dithiolane",
                        "benzothiophene",
                        "oxathiolane",
                        "dioxathiolane",
                        "thiazolidine",
                        "oxazolidine",
                        "isoxazole",
                        "isothiazole",
                        "oxadiazole",
                        "thiadiazole",
                        "benzoxazole",
                        "benzothiazole",
                        "benzimidazole",
                        "pteridin",
                        "phenothiazine",
                        "phenoxazine",
                        "dibenzofuran",
                        "dibenzothiophene",
                        "xanthene",
                        "thioxanthene",
                        "pyrroline",
                        "pyrrolidone",
                        "imidazolidine",
                        "porphyrin",
                        "indazole",
                        "benzotriazole",
                    ]

                    # Check if any heterocycle is formed in this reaction
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product):
                            # Check if the ring was not present in the reactants
                            reactants = rsmi.split(">")[0].split(".")
                            ring_in_reactants = any(
                                checker.check_ring(ring, reactant) for reactant in reactants
                            )

                            if not ring_in_reactants:
                                print(f"Heterocycle formation detected: {ring} at depth {depth}")
                                heterocycle_depths.append(depth)
                except Exception as e:
                    print(f"Error in heterocycle detection: {e}")

            # Recursively process children
            for child in node.get("children", []):
                dfs_heterocycle(child, depth + 1)

        # Start DFS traversal
        dfs_heterocycle(route)

        # Check if heterocycles are formed in early stages (higher depth values)
        if heterocycle_depths:
            max_depth = max(heterocycle_depths)
            total_depth = get_max_depth(route)
            # Consider it early stage if in the first half of the synthesis
            is_early_stage = max_depth > total_depth / 2
            print(
                f"Heterocycle formation depths: {heterocycle_depths}, max depth: {max_depth}, total depth: {total_depth}"
            )
            return is_early_stage
        return False

    def late_stage_ether_formation_strategy(route):
        """Check if ether formation occurs in late stages of synthesis."""
        ether_depths = []

        def dfs_ether(node, depth=0):
            if node["type"] == "reaction":
                try:
                    rsmi = node["metadata"]["rsmi"]

                    # Check for ether formation reactions
                    ether_reactions = [
                        "Williamson Ether Synthesis",
                        "Williamson Ether Synthesis (intra to epoxy)",
                        "Mitsunobu aryl ether",
                        "Mitsunobu aryl ether (intramolecular)",
                        "{Williamson ether}",
                        "Chan-Lam etherification",
                        "Alcohol to ether",
                        "Ullmann-Goldberg Substitution aryl alcohol",
                    ]

                    for rxn_type in ether_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Ether formation detected: {rxn_type} at depth {depth}")
                            ether_depths.append(depth)
                            break

                    # Also check for ether functional group formation if not detected by reaction type
                    if not ether_depths or ether_depths[-1] != depth:
                        product = rsmi.split(">")[-1]
                        reactants = rsmi.split(">")[0].split(".")

                        # Check if number of ethers increased
                        ether_in_product = len(checker.get_fg_atom_indices("Ether", product))
                        ether_in_reactants = sum(
                            len(checker.get_fg_atom_indices("Ether", r)) for r in reactants if r
                        )

                        if ether_in_product > ether_in_reactants:
                            print(
                                f"Ether formation detected by functional group count at depth {depth}"
                            )
                            ether_depths.append(depth)

                        # Additional check for specific ether-forming reactions not covered above
                        if checker.check_reaction(
                            "Ullmann condensation", rsmi
                        ) and checker.check_fg("Ether", product):
                            print(
                                f"Ether formation detected via Ullmann condensation at depth {depth}"
                            )
                            ether_depths.append(depth)
                except Exception as e:
                    print(f"Error in ether formation detection: {e}")

            # Recursively process children
            for child in node.get("children", []):
                dfs_ether(child, depth + 1)

        # Start DFS traversal
        dfs_ether(route)

        # Check if ether formation occurs in late stages (lower depth values)
        if ether_depths:
            min_depth = min(ether_depths)
            total_depth = get_max_depth(route)
            # Consider it late stage if in the first third of the synthesis
            is_late_stage = min_depth <= total_depth / 3
            print(
                f"Ether formation depths: {ether_depths}, min depth: {min_depth}, total depth: {total_depth}"
            )
            return is_late_stage
        return False

    def get_max_depth(node, depth=0):
        """Helper function to get the maximum depth of the synthesis tree."""
        if not node.get("children", []):
            return depth

        max_child_depth = 0
        for child in node.get("children", []):
            child_depth = get_max_depth(child, depth + 1)
            max_child_depth = max(max_child_depth, child_depth)

        return max_child_depth

    # Use the individual strategy functions
    is_heterocycle_first = heterocycle_first_strategy(route)
    is_late_stage_ether = late_stage_ether_formation_strategy(route)

    # Combined strategy
    combined_strategy = is_heterocycle_first and is_late_stage_ether

    print(f"Heterocycle-first strategy detected: {is_heterocycle_first}")
    print(f"Late-stage ether formation strategy detected: {is_late_stage_ether}")
    print(f"Heterocycle-first with late-stage ether strategy detected: {combined_strategy}")

    return combined_strategy
