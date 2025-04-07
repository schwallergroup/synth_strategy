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
    This function detects if the synthesis follows a linear strategy of
    functionalizing a piperidine scaffold at multiple positions.
    """
    # Track piperidine-containing molecules at different depths
    piperidine_steps = []
    functionalization_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains a piperidine ring
            if checker.check_ring("piperidine", mol_smiles):
                print(f"Found piperidine at depth {depth}: {mol_smiles}")
                piperidine_steps.append((depth, mol_smiles))

            # Continue traversal
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

        elif node["type"] == "reaction":
            # Extract reaction information
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                # Continue traversal if no reaction SMILES
                for child in node.get("children", []):
                    dfs_traverse(child, depth)
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains piperidine
            if checker.check_ring("piperidine", product):
                # Check if any reactant also contains piperidine
                piperidine_reactant = None
                for reactant in reactants:
                    if checker.check_ring("piperidine", reactant):
                        piperidine_reactant = reactant
                        break

                # If no piperidine in reactants but in product, it's a piperidine formation
                if not piperidine_reactant:
                    print(f"Found piperidine formation at depth {depth}")
                    print(f"Reaction: {rsmi}")
                    functionalization_events.append((depth, rsmi))

                # If we found a piperidine in both product and reactant, check for functionalization
                elif piperidine_reactant:
                    # Check for common functionalization reactions
                    if (
                        checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Acylation of secondary amines", rsmi)
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Boc amine protection", rsmi)
                        or checker.check_reaction("Boc amine deprotection", rsmi)
                        or checker.check_reaction("N-arylation", rsmi)
                    ):

                        print(f"Found piperidine functionalization reaction at depth {depth}")
                        print(f"Reaction: {rsmi}")
                        functionalization_events.append((depth, rsmi))

                    # Check for other common functional group additions
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mol = Chem.MolFromSmiles(piperidine_reactant)

                    if product_mol and reactant_mol:
                        # Check for new functional groups in the product
                        new_functionalization = False
                        for fg in [
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Ester",
                            "Primary amide",
                            "Secondary amide",
                            "Tertiary amide",
                            "Nitrile",
                            "Ketone",
                            "Aldehyde",
                            "Carboxylic acid",
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                            "Aromatic halide",
                            "Sulfonamide",
                            "Urea",
                            "Carbamate",
                            "Primary amine",
                            "Secondary amine",
                            "Tertiary amine",
                        ]:
                            if checker.check_fg(fg, product) and not checker.check_fg(
                                fg, piperidine_reactant
                            ):
                                print(f"Added {fg} to piperidine at depth {depth}")
                                new_functionalization = True

                        if new_functionalization:
                            functionalization_events.append((depth, rsmi))

                        # If the piperidine structure changed but we didn't catch it with specific checks
                        elif piperidine_reactant != product:
                            print(f"Piperidine structure changed at depth {depth}")
                            print(f"From: {piperidine_reactant}")
                            print(f"To: {product}")
                            functionalization_events.append((depth, rsmi))

            # Check if any reactant contains piperidine but product doesn't
            elif any(checker.check_ring("piperidine", reactant) for reactant in reactants):
                print(f"Piperidine transformation at depth {depth}")
                print(f"Reaction: {rsmi}")
                functionalization_events.append((depth, rsmi))

            # Continue traversal
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze results
    print(f"Found piperidine at {len(piperidine_steps)} steps")
    print(f"Found {len(functionalization_events)} functionalization events")

    # Check if we have a linear piperidine functionalization strategy
    # Need at least 3 steps with piperidine and at least 2 functionalization events
    if len(piperidine_steps) >= 3 and len(functionalization_events) >= 2:
        # Check if the functionalizations occur at different depths
        functionalization_depths = set(depth for depth, _ in functionalization_events)
        if len(functionalization_depths) >= 2:
            print("Found linear piperidine functionalization strategy")
            return True

    return False
