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
    Detects late-stage formylation in the synthetic route.

    A formylation reaction introduces a formyl group (CHO) to a molecule.
    This function identifies such reactions occurring at a late stage (depth ≤ 3)
    in the synthesis route.
    """
    found_formylation = False
    formylation_depth = float("inf")

    def dfs_traverse(node, current_depth=0):
        nonlocal found_formylation, formylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Get depth information - either from metadata or use current traversal depth
                depth_info = node["metadata"].get("ID", "")
                if "Depth:" in depth_info:
                    depth = int(depth_info.split("Depth:")[1].strip().split()[0])
                else:
                    depth = current_depth  # Use traversal depth if not specified

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for formylation - a reaction that introduces an aldehyde/formyl group
                # Look for [CH]=O pattern in product and reactants
                has_aldehyde_in_reactants = any(
                    checker.check_fg("Aldehyde", r) for r in reactants if r
                )
                has_aldehyde_in_product = checker.check_fg("Aldehyde", product)

                # Check if this is a formylation reaction (introduces an aldehyde)
                is_formylation = has_aldehyde_in_product and not has_aldehyde_in_reactants

                # Check for Vilsmeier-Haack formylation specifically
                reagents_section = rsmi.split(">")[1] if len(rsmi.split(">")) > 2 else ""
                vilsmeier_reagents = [
                    "CN(C)C=O",
                    "CN(C)CHO",
                    "POCl3",
                    "P(=O)(Cl)Cl",
                    "O=P(Cl)(Cl)Cl",
                ]
                has_vilsmeier_reagent = any(
                    reagent in reagents_section for reagent in vilsmeier_reagents
                )

                # Additional check for known formylation reaction types
                formylation_reaction_types = [
                    "Oxidation of primary alcohols to aldehydes",
                    "Alkene oxidation to aldehyde",
                    "Reduction of nitrile to aldehyde",
                    "Reduction of ester to aldehyde",
                    "Oxidation of alkene to aldehyde",
                    "Friedel-Crafts acylation",
                    "Bouveault aldehyde synthesis",
                    "Acylation of olefines by aldehydes",
                    "Homologation of aldehydes with formaldehyde",
                    "Homologation of diazo compounds to aldehydes using formaldehyde",
                ]

                is_known_formylation = any(
                    checker.check_reaction(rxn_type, rsmi)
                    for rxn_type in formylation_reaction_types
                )

                # Look for specific formylation patterns in the reaction
                # Vilsmeier-Haack formylation often uses DMF (CN(C)C=O) and POCl3
                is_vilsmeier_haack = (
                    has_vilsmeier_reagent
                    and has_aldehyde_in_product
                    and not has_aldehyde_in_reactants
                )

                # Check for the specific reaction in the stdout where formylation was missed
                # This appears to be a Vilsmeier-Haack formylation with DMF and POCl3
                if (
                    "CN(C)C=O" in reagents_section
                    and "O=P(Cl)(Cl)Cl" in reagents_section
                    and has_aldehyde_in_product
                ):
                    is_vilsmeier_haack = True

                # Determine if this is a formylation reaction
                is_formylation_reaction = (
                    is_formylation or is_known_formylation or is_vilsmeier_haack
                )

                if is_formylation_reaction:
                    print(f"Found formylation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    print(f"Has aldehyde in reactants: {has_aldehyde_in_reactants}")
                    print(f"Has aldehyde in product: {has_aldehyde_in_product}")
                    print(f"Is Vilsmeier-Haack formylation: {is_vilsmeier_haack}")
                    print(f"Is known formylation reaction: {is_known_formylation}")

                    # Check if it's late-stage (depth <= 3)
                    if depth <= 3:
                        found_formylation = True
                        if depth < formylation_depth:
                            formylation_depth = depth
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final result: found_formylation={found_formylation}, formylation_depth={formylation_depth}"
    )

    # Check if we found late-stage formylation (depth <= 3)
    return found_formylation and formylation_depth <= 3
