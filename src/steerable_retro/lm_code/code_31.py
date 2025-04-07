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
    This function detects if the synthesis involves late-stage O-alkylation
    (defined as O-alkylation in the final step or depth 0 or 1).
    """
    late_stage_o_alkylation_detected = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_o_alkylation_detected

        # Set depth for the current node
        node["depth"] = current_depth

        print(f"Traversing node at depth {current_depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            depth = node.get("depth", None)
            if depth <= 1:  # Late stage (depth 0 or 1)
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for various O-alkylation reactions
                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: Williamson Ether Synthesis at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: Williamson Ether Synthesis (intra to epoxy) at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                ):
                    print(
                        f"Detected late-stage O-alkylation: O-alkylation of carboxylic acids with diazo compounds at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("O-alkylation of amides with diazo compounds", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: O-alkylation of amides with diazo compounds at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Mitsunobu aryl ether", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: Mitsunobu aryl ether at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("{Williamson ether}", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: {{Williamson ether}} at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Mitsunobu_phenole", rsmi):
                    print(f"Detected late-stage O-alkylation: Mitsunobu_phenole at depth {depth}")
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Chan-Lam etherification", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: Chan-Lam etherification at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("O-methylation", rsmi):
                    print(f"Detected late-stage O-alkylation: O-methylation at depth {depth}")
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Methylation of OH with DMS", rsmi):
                    print(
                        f"Detected late-stage O-alkylation: Methylation of OH with DMS at depth {depth}"
                    )
                    late_stage_o_alkylation_detected = True
                elif checker.check_reaction("Alcohol to ether", rsmi):
                    print(f"Detected late-stage O-alkylation: Alcohol to ether at depth {depth}")
                    late_stage_o_alkylation_detected = True
                else:
                    # Check for general O-alkylation pattern if specific reaction types aren't detected
                    has_alcohol_reactant = False
                    has_alkylating_agent = False

                    for reactant in reactants:
                        # Check if reactant contains alcohol or phenol
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Phenol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                        ):
                            has_alcohol_reactant = True
                            print(f"Found alcohol/phenol in reactant: {reactant}")

                        # Check if reactant is an alkylating agent
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Triflate", reactant)
                            or checker.check_fg("Mesylate", reactant)
                            or checker.check_fg("Tosylate", reactant)
                        ):
                            has_alkylating_agent = True
                            print(f"Found alkylating agent in reactant: {reactant}")

                    # Check if product contains ether
                    has_ether_product = checker.check_fg("Ether", product)
                    if has_ether_product:
                        print(f"Found ether in product: {product}")

                    # Verify O-alkylation pattern
                    if has_alcohol_reactant and has_alkylating_agent and has_ether_product:
                        print(f"Detected late-stage O-alkylation: General pattern at depth {depth}")
                        late_stage_o_alkylation_detected = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_stage_o_alkylation_detected}")
    return late_stage_o_alkylation_detected
