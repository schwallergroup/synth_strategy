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
    This function detects multiple N-alkylation steps on a piperazine scaffold.
    It looks for at least two separate reactions where a nitrogen in a piperazine
    ring is alkylated.
    """
    # Track alkylation reactions on piperazine scaffolds
    piperazine_alkylations = {}
    target_mol = None

    def dfs_traverse(node, depth=0, path=None):
        nonlocal target_mol

        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Depth {depth}, Examining molecule: {mol_smiles[:30]}...")

            # If this molecule contains piperazine, it could be our target
            if checker.check_ring("piperazine", mol_smiles):
                print(f"Found piperazine in molecule at depth {depth}")

                # If we don't have a target yet, or this is a later stage (lower depth)
                if target_mol is None or depth < target_mol[1]:
                    target_mol = (mol_smiles, depth)
                    print(f"Set target molecule to: {mol_smiles[:30]}... at depth {depth}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth}, Examining reaction: {rsmi[:50]}...")

                # Check if product contains piperazine
                if checker.check_ring("piperazine", product):
                    print(f"Product contains piperazine")

                    # Find reactant with piperazine
                    piperazine_reactant = None
                    for reactant in reactants:
                        if checker.check_ring("piperazine", reactant):
                            piperazine_reactant = reactant
                            print(f"Found piperazine in reactant: {reactant[:30]}...")
                            break

                    if piperazine_reactant:
                        # Check if this is an N-alkylation reaction
                        n_alkylation_reactions = [
                            "N-alkylation of primary amines with alkyl halides",
                            "N-alkylation of secondary amines with alkyl halides",
                            "Methylation with MeI_primary",
                            "Methylation with MeI_secondary",
                            "Methylation with MeI_tertiary",
                            "Methylation",
                            "DMS Amine methylation",
                            "Eschweiler-Clarke Primary Amine Methylation",
                            "Eschweiler-Clarke Secondary Amine Methylation",
                            "Reductive methylation of primary amine with formaldehyde",
                            "N-methylation",
                            "Reductive amination with aldehyde",
                            "Reductive amination with ketone",
                            "Reductive amination with alcohol",
                            "Alkylation of amines",
                        ]

                        for reaction_type in n_alkylation_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(f"Detected N-alkylation reaction: {reaction_type}")

                                # Store this alkylation with the product as key
                                if product not in piperazine_alkylations:
                                    piperazine_alkylations[product] = []
                                piperazine_alkylations[product].append((reaction_type, depth))
                                print(
                                    f"Added alkylation to product {product[:30]}... at depth {depth}"
                                )
                                break

                        # If no specific reaction type was found, check if there's evidence of N-alkylation
                        # by comparing secondary/tertiary amines in reactant vs product
                        if product not in piperazine_alkylations:
                            reactant_has_secondary = checker.check_fg(
                                "Secondary amine", piperazine_reactant
                            )
                            reactant_has_tertiary = checker.check_fg(
                                "Tertiary amine", piperazine_reactant
                            )
                            product_has_tertiary = checker.check_fg("Tertiary amine", product)

                            if (reactant_has_secondary and product_has_tertiary) or (
                                reactant_has_tertiary
                                and not checker.check_fg("Secondary amine", product)
                            ):
                                print(f"Detected N-alkylation by functional group change")
                                if product not in piperazine_alkylations:
                                    piperazine_alkylations[product] = []
                                piperazine_alkylations[product].append(
                                    ("N-alkylation (detected by FG change)", depth)
                                )
                                print(
                                    f"Added alkylation to product {product[:30]}... at depth {depth}"
                                )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Count unique alkylation reactions on the target scaffold
    alkylation_count = 0
    if target_mol:
        target_smiles, _ = target_mol
        print(f"Target molecule with piperazine: {target_smiles[:30]}...")

        # Count alkylations in the synthetic path to the target
        for product, alkylations in piperazine_alkylations.items():
            # Check if this product is in the path to our target
            if checker.check_ring("piperazine", product):
                for alkylation_info in alkylations:
                    alkylation_count += 1
                    print(
                        f"Counting alkylation: {alkylation_info[0]} at depth {alkylation_info[1]}"
                    )

    print(f"Total N-alkylation count on piperazine scaffold: {alkylation_count}")
    return alkylation_count >= 2
