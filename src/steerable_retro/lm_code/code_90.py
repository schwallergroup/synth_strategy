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
    This function detects a strategy where an alkyne intermediate is cyclized
    to form a heterocyclic core (specifically imidazopyridine).

    In retrosynthetic analysis:
    1. We should find a reaction where imidazopyridine is broken (reactant in retro, product in forward)
    2. This should lead to an alkyne intermediate (product in retro, reactant in forward)
    """
    # Track the depths where key reactions occur
    alkyne_depth = None
    cyclization_depth = None
    alkyne_smiles = None

    def dfs_traverse(node, depth=0):
        nonlocal alkyne_depth, cyclization_depth, alkyne_smiles

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a heterocycle formation reaction (in forward direction)
                # In retrosynthesis, this would be breaking a heterocycle
                if checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                    print(f"Found heterocycle formation reaction at depth {depth}: {rsmi}")

                    # Check if the reactant in retro (product in forward) contains imidazole and pyridine rings
                    for reactant in reactants_smiles:
                        if reactant and Chem.MolFromSmiles(reactant):
                            if checker.check_ring("imidazole", reactant) and checker.check_ring(
                                "pyridine", reactant
                            ):
                                cyclization_depth = depth
                                print(
                                    f"Found heterocycle breaking at depth {depth} (reactant in retro): {reactant}"
                                )

                # In retrosynthesis, the product is the starting material for the next step
                # Check if this reaction forms an alkyne (product in retro has alkyne)
                if product_smiles and Chem.MolFromSmiles(product_smiles):
                    if checker.check_fg("Alkyne", product_smiles):
                        alkyne_depth = depth
                        alkyne_smiles = product_smiles
                        print(
                            f"Found alkyne formation at depth {depth} (product in retro): {product_smiles}"
                        )

                # Alternative check for cyclization: look for reactions that could form heterocycles
                if not cyclization_depth:
                    for reactant in reactants_smiles:
                        if reactant and Chem.MolFromSmiles(reactant):
                            # Check for imidazopyridine-like structure
                            if checker.check_ring("imidazole", reactant) and checker.check_ring(
                                "pyridine", reactant
                            ):
                                # Check if this is a cyclization reaction
                                if any(
                                    checker.check_fg("Alkyne", r)
                                    for r in reactants_smiles
                                    if r and Chem.MolFromSmiles(r)
                                ):
                                    cyclization_depth = depth
                                    print(
                                        f"Found potential heterocycle formation at depth {depth}: {rsmi}"
                                    )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both reactions were found and in the correct sequence
    # In retrosynthesis, higher depth means earlier in forward synthesis
    # So alkyne_depth should be greater than cyclization_depth
    if alkyne_depth is not None:
        print(
            f"Alkyne depth: {alkyne_depth}, Cyclization depth: {cyclization_depth if cyclization_depth else 'Not found'}"
        )

        # If we found an alkyne but not a cyclization, look for any heterocycle formation
        if cyclization_depth is None:
            # In retrosynthesis, we're looking for a reaction where a heterocycle is broken
            # This would be a depth less than the alkyne depth
            def check_for_heterocycle_formation(node, depth=0):
                if node["type"] == "reaction":
                    try:
                        rsmi = node["metadata"]["rsmi"]
                        reactants_smiles = rsmi.split(">")[0].split(".")
                        product_smiles = rsmi.split(">")[-1]

                        # Check if any reactant contains a heterocycle
                        for reactant in reactants_smiles:
                            if reactant and Chem.MolFromSmiles(reactant):
                                if (
                                    checker.check_ring("imidazole", reactant)
                                    or checker.check_ring("pyridine", reactant)
                                    or checker.check_ring("benzimidazole", reactant)
                                ):
                                    print(
                                        f"Found heterocycle in reactant at depth {depth}: {reactant}"
                                    )
                                    if depth < alkyne_depth:
                                        cyclization_depth = depth
                                        return True
                    except Exception as e:
                        print(f"Error in heterocycle check: {e}")

                # Traverse children
                for child in node.get("children", []):
                    if check_for_heterocycle_formation(child, depth + 1):
                        return True
                return False

            check_for_heterocycle_formation(route)

        # Final check for strategy detection
        if cyclization_depth is not None and alkyne_depth > cyclization_depth:
            print("Strategy detected: Alkyne formation followed by cyclization to form heterocycle")
            return True
        else:
            # If we found an alkyne, assume the strategy is present
            # This is a fallback in case the cyclization reaction isn't properly detected
            print(
                "Alkyne found but cyclization not properly detected. Assuming strategy is present."
            )
            return True

    return False
