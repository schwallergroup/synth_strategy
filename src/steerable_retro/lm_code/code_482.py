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
    Detects a strategy involving sequential functionalization of a heterocyclic scaffold
    at multiple positions, maintaining the core structure throughout.
    """
    # Track key features
    functionalization_steps = 0
    heterocycle_type = None
    different_positions_modified = set()

    # List of heterocycles to check
    heterocycle_types = [
        "furan",
        "pyran",
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
        "indole",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of functional groups to check for modifications
    functional_groups = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aniline",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Phenol",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Nitro group",
        "Aldehyde",
        "Ketone",
        "Ether",
        "Thiol",
        "Sulfide",
        "Sulfoxide",
        "Sulfone",
    ]

    # List of common heterocycle functionalization reactions
    functionalization_reactions = [
        "Aromatic bromination",
        "Aromatic chlorination",
        "Aromatic fluorination",
        "Aromatic iodination",
        "Suzuki coupling",
        "Buchwald-Hartwig",
        "N-arylation",
        "Heck terminal vinyl",
        "Sonogashira",
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Minisci",
        "Directed ortho metalation of arenes",
        "Nucleophilic substitution",
    ]

    def identify_heterocycle(smiles):
        """Identify which heterocycle is present in the molecule"""
        for hc_type in heterocycle_types:
            if checker.check_ring(hc_type, smiles):
                print(f"Found heterocycle: {hc_type} in {smiles}")
                return hc_type
        return None

    def is_functionalization_reaction(rxn_smiles):
        """Check if the reaction is a common heterocycle functionalization reaction"""
        for rxn_type in functionalization_reactions:
            if checker.check_reaction(rxn_type, rxn_smiles):
                print(f"Detected functionalization reaction: {rxn_type}")
                return True
        return False

    def detect_functionalization(reactants_smiles, product_smiles, rxn_smiles, hc_type):
        """
        Detect if a functionalization occurred on the heterocycle by comparing
        functional groups in reactants and products
        """
        if not hc_type:
            return set()

        try:
            # First check if this is a known functionalization reaction
            if is_functionalization_reaction(rxn_smiles):
                # If it's a known functionalization reaction and the heterocycle is preserved,
                # we'll consider it a valid functionalization step
                if checker.check_ring(hc_type, reactants_smiles) and checker.check_ring(
                    hc_type, product_smiles
                ):
                    # For simplicity, we'll just return a placeholder position
                    # This is a simplification since we can't reliably determine the exact position
                    return {1}

            # Check for functional group changes
            reactant_fgs = set()
            product_fgs = set()

            for fg in functional_groups:
                if checker.check_fg(fg, reactants_smiles):
                    reactant_fgs.add(fg)
                if checker.check_fg(fg, product_smiles):
                    product_fgs.add(fg)

            # If there are differences in functional groups and the heterocycle is preserved,
            # consider it a functionalization step
            if (
                reactant_fgs != product_fgs
                and checker.check_ring(hc_type, reactants_smiles)
                and checker.check_ring(hc_type, product_smiles)
            ):
                added_fgs = product_fgs - reactant_fgs
                removed_fgs = reactant_fgs - product_fgs

                if added_fgs or removed_fgs:
                    print(f"Functional group changes: Added {added_fgs}, Removed {removed_fgs}")
                    # For simplicity, we'll just return a unique position for each change
                    # This is a simplification since we can't reliably determine the exact position
                    return {hash(fg) % 100 for fg in (added_fgs | removed_fgs)}

            return set()
        except Exception as e:
            print(f"Error in detect_functionalization: {e}")
            return set()

    def dfs_traverse(node, depth=0):
        nonlocal functionalization_steps, heterocycle_type, different_positions_modified

        if node["type"] == "mol":
            # Check for heterocycle in molecules
            mol_smiles = node["smiles"]
            current_heterocycle = identify_heterocycle(mol_smiles)

            # Initialize heterocycle_type if not set
            if heterocycle_type is None and current_heterocycle:
                heterocycle_type = current_heterocycle
                print(f"Initial heterocycle type set to: {heterocycle_type}")

        elif node["type"] == "reaction" and heterocycle_type:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Verify the heterocycle is present in both reactants and product
                if checker.check_ring(heterocycle_type, reactants) and checker.check_ring(
                    heterocycle_type, product
                ):
                    # Check if this reaction modifies the heterocycle
                    modified_positions = detect_functionalization(
                        reactants, product, rsmi, heterocycle_type
                    )

                    if modified_positions:
                        functionalization_steps += 1
                        different_positions_modified.update(modified_positions)
                        print(
                            f"Found functionalization step modifying positions: {modified_positions}"
                        )
                        print(f"Reaction SMILES: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        functionalization_steps >= 2
        and heterocycle_type is not None
        and len(different_positions_modified) >= 2
    )

    print(f"Heterocycle sequential functionalization strategy detected: {strategy_present}")
    print(f"Heterocycle type: {heterocycle_type}")
    print(
        f"Functionalization steps: {functionalization_steps}, Different positions: {len(different_positions_modified)}"
    )

    return strategy_present
