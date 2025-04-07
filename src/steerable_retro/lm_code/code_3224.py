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
    This function detects a linear synthetic strategy for constructing heterocycles
    with no convergent steps and a final cyclization.
    """
    # Initialize tracking variables
    reaction_count = 0
    has_final_cyclization = False
    max_reactants_per_step = 0
    heterocycle_formed = False
    final_product_smiles = ""

    # List of common heterocycles to check
    heterocycle_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
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
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # List of cyclization reaction types
    cyclization_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "imidazole",
        "Intramolecular amination (heterocycle formation)",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
    ]

    # Get the final product SMILES from the root node
    if route["type"] == "mol":
        final_product_smiles = route["smiles"]
        print(f"Final product: {final_product_smiles}")

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, has_final_cyclization, max_reactants_per_step, heterocycle_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                print(f"Invalid reaction SMILES format at depth {depth}: {rsmi}")
                return

            reactants_smiles = [r for r in parts[0].split(".") if r]
            product_smiles = parts[2]

            print(f"Processing reaction at depth {depth}: {rsmi}")

            # Count this reaction
            reaction_count += 1

            # Track maximum number of reactants in any step
            max_reactants_per_step = max(max_reactants_per_step, len(reactants_smiles))
            print(f"Reactants in this step: {len(reactants_smiles)}")

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if not product_mol or not all(reactant_mols):
                    print("Invalid molecules detected, skipping reaction")
                    return

                # Check for ring formation in final step (depth 1 in retrosynthetic tree)
                if depth == 1:
                    print(f"Analyzing final synthetic step at depth {depth}")

                    # Check if this is a known cyclization reaction
                    for rxn_type in cyclization_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected cyclization reaction: {rxn_type}")
                            has_final_cyclization = True
                            break

                    # Count rings in product and reactants
                    product_ring_count = len(Chem.GetSSSR(product_mol))
                    reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                    max_reactant_ring_count = max(reactant_ring_counts, default=0)

                    print(
                        f"Final step - Product rings: {product_ring_count}, Max reactant rings: {max_reactant_ring_count}"
                    )

                    # Check if a new ring was formed
                    if product_ring_count > max_reactant_ring_count:
                        print("New ring formed in final step")

                        # Check if the product contains a heterocycle
                        for heterocycle in heterocycle_types:
                            if checker.check_ring(heterocycle, product_smiles):
                                # Check if this heterocycle wasn't in any reactant
                                heterocycle_in_reactants = any(
                                    checker.check_ring(heterocycle, r) for r in reactants_smiles
                                )

                                if not heterocycle_in_reactants:
                                    print(f"Heterocycle {heterocycle} formed in final step")
                                    has_final_cyclization = True
                                    heterocycle_formed = True
                                    break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root node
    dfs_traverse(route)

    # Strategy is present if we have multiple steps, final cyclization, and no convergent steps
    # (max 2 reactants per step indicates linear synthesis)
    strategy_present = (
        reaction_count >= 3
        and has_final_cyclization
        and max_reactants_per_step <= 2
        and heterocycle_formed
    )

    print(f"Linear heterocycle construction strategy detected: {strategy_present}")
    print(f"Reaction count: {reaction_count}, Max reactants per step: {max_reactants_per_step}")
    print(f"Final cyclization: {has_final_cyclization}, Heterocycle formed: {heterocycle_formed}")

    return strategy_present
