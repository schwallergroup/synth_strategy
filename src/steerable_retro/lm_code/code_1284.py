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
    Detects N-alkylation as a key bond-forming step in the synthesis.
    """
    # Track if we found N-alkylation
    found_n_alkylation = False
    target_mol = route["smiles"]
    print(f"Target molecule: {target_mol}")

    def dfs_traverse(node, depth=0):
        nonlocal found_n_alkylation

        if node["type"] == "reaction" and depth <= 3:  # Focus on late-stage reactions
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES at depth {depth}")
                    return

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for N-alkylation reaction types directly
                n_alkylation_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "Methylation",
                    "N-methylation",
                    "Eschweiler-Clarke Primary Amine Methylation",
                    "Eschweiler-Clarke Secondary Amine Methylation",
                    "Reductive methylation of primary amine with formaldehyde",
                ]

                for rxn_type in n_alkylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found N-alkylation reaction ({rxn_type}) at depth {depth}: {rsmi}")
                        found_n_alkylation = True
                        return

                # If specific reaction check fails, try checking reactants and products
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if any reactant is an amine (primary or secondary)
                has_amine = False
                for reactant in reactants_str.split("."):
                    if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                        "Secondary amine", reactant
                    ):
                        has_amine = True
                        print(f"Found amine in reactant: {reactant}")
                        break

                # Check if any reactant is an alkyl halide or carbonyl (for reductive amination)
                has_alkylating_agent = False
                for reactant in reactants_str.split("."):
                    if (
                        checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                        or checker.check_fg("Tertiary halide", reactant)
                        or checker.check_fg("Aldehyde", reactant)
                        or checker.check_fg("Ketone", reactant)
                    ):
                        has_alkylating_agent = True
                        print(f"Found alkylating agent in reactant: {reactant}")
                        break

                # Check if product has a more substituted amine
                has_more_substituted_amine = False
                if checker.check_fg("Secondary amine", product_str) or checker.check_fg(
                    "Tertiary amine", product_str
                ):
                    has_more_substituted_amine = True
                    print(f"Found more substituted amine in product: {product_str}")

                # If we have an amine, an alkylating agent, and the product has a more substituted amine,
                # it's likely an N-alkylation
                if has_amine and has_alkylating_agent and has_more_substituted_amine:
                    print(
                        f"Found N-alkylation by functional group analysis at depth {depth}: {rsmi}"
                    )
                    found_n_alkylation = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"N-alkylation found: {found_n_alkylation}")

    return found_n_alkylation
