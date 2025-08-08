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
    Detects a late-stage convergent synthesis strategy where two halogenated
    heterocyclic fragments are coupled in the final or penultimate step.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track halogenated heterocycles
    halogenated_heterocycles = []
    # Track final coupling
    final_coupling = None

    # List of heterocyclic rings to check
    heterocycle_rings = [
        "thiophene",
        "pyridine",
        "triazole",
        "furan",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "indole",
        "benzothiophene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, halogenated_heterocycles, final_coupling

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this is a halogenation reaction
                is_halogenation = (
                    checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                    or checker.check_reaction("Bromination", rsmi)
                    or checker.check_reaction("Iodination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Chlorination", rsmi)
                )

                # Check if product contains a heterocycle and a halogen regardless of reaction type
                has_heterocycle = False
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, product_part):
                        has_heterocycle = True
                        print(f"Found heterocycle {ring} in product at depth {depth}")
                        break

                has_halogen = "Br" in product_part or "I" in product_part or "Cl" in product_part

                if has_heterocycle and has_halogen:
                    print(f"Found halogenated heterocycle at depth {depth}: {product_part}")
                    halogenated_heterocycles.append((depth, product_part))

                # Check if this is a coupling reaction
                is_coupling = (
                    checker.check_reaction("Suzuki", rsmi)
                    or checker.check_reaction("Stille", rsmi)
                    or checker.check_reaction("Negishi", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                    or checker.check_reaction("Sonogashira", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution", rsmi)
                )

                # If not a known coupling, check for any cross-coupling-like reaction
                if not is_coupling and len(reactants) >= 2:
                    has_halogen = any("Br" in r or "I" in r or "Cl" in r for r in reactants)
                    has_heterocycle_reactants = sum(
                        1
                        for r in reactants
                        if any(checker.check_ring(ring, r) for ring in heterocycle_rings)
                    )

                    if has_halogen and has_heterocycle_reactants >= 1:
                        print(f"Detected potential custom coupling reaction at depth {depth}")
                        is_coupling = True

                if is_coupling and depth <= 1:
                    # Check if reactants have heterocycles
                    heterocycle_count = 0
                    halogen_count = 0

                    for reactant in reactants:
                        has_heterocycle = False
                        for ring in heterocycle_rings:
                            if checker.check_ring(ring, reactant):
                                has_heterocycle = True
                                print(
                                    f"Found heterocycle {ring} in coupling reactant at depth {depth}: {reactant}"
                                )
                                heterocycle_count += 1
                                break

                        has_halogen = "Br" in reactant or "I" in reactant or "Cl" in reactant
                        if has_halogen:
                            halogen_count += 1
                            print(
                                f"Found halogen in coupling reactant at depth {depth}: {reactant}"
                            )

                    if heterocycle_count >= 1 and halogen_count >= 1:
                        print(f"Found potential coupling at depth {depth}: {rsmi}")
                        final_coupling = (depth, rsmi, reactants)
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both halogenated heterocycles and a final coupling
    if halogenated_heterocycles and final_coupling is not None:
        coupling_depth, coupling_rsmi, coupling_reactants = final_coupling

        # Ensure coupling is late-stage (depth 0 or 1)
        if coupling_depth <= 1:
            # Count heterocycles in the coupling
            heterocycle_count_in_coupling = 0
            halogenated_heterocycle_in_coupling = False

            for reactant in coupling_reactants:
                # Check if this reactant is a heterocycle
                is_heterocycle = False
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, reactant):
                        is_heterocycle = True
                        heterocycle_count_in_coupling += 1
                        print(f"Confirmed heterocycle {ring} in coupling reactant")
                        break

                # Check if this reactant has a halogen
                has_halogen = "Br" in reactant or "I" in reactant or "Cl" in reactant

                # If it's a halogenated heterocycle
                if is_heterocycle and has_halogen:
                    halogenated_heterocycle_in_coupling = True
                    print(f"Found halogenated heterocycle in coupling: {reactant}")

            # For convergent synthesis, we need at least one halogenated heterocycle in a late-stage coupling
            if (
                coupling_depth <= 1
                and heterocycle_count_in_coupling >= 1
                and halogenated_heterocycle_in_coupling
            ):
                found_pattern = True
                print("Confirmed late-stage convergent halogenated heterocycle coupling strategy")

    return found_pattern
