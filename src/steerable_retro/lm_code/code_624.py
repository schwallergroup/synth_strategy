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
    This function detects a strategy involving indole ring formation followed by
    N-sulfonylation and late-stage oxidation of a benzylic position.
    """
    # Initialize flags for each key transformation
    indole_formation_detected = False
    n_sulfonylation_detected = False
    benzylic_modification_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal indole_formation_detected, n_sulfonylation_detected, benzylic_modification_detected

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found")
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Reaction at depth {depth}: {rsmi}")
            print(f"Product: {product_smiles}")
            print(f"Reactants: {reactants_smiles}")

            try:
                # Check for indole formation (can occur at any depth >= 2)
                if not indole_formation_detected and depth >= 2:
                    # Check if product contains indole but reactants don't
                    if checker.check_ring("indole", product_smiles) and not any(
                        checker.check_ring("indole", r) for r in reactants_smiles
                    ):
                        print(f"Detected indole formation at depth {depth}")
                        indole_formation_detected = True
                    # Check for Fischer indole synthesis
                    elif checker.check_reaction("Fischer indole", rsmi):
                        print(f"Detected Fischer indole synthesis at depth {depth}")
                        indole_formation_detected = True
                    else:
                        print(f"No indole formation detected at depth {depth}")
                        print(f"Product has indole: {checker.check_ring('indole', product_smiles)}")
                        reactants_with_indole = [
                            r for r in reactants_smiles if checker.check_ring("indole", r)
                        ]
                        print(f"Reactants with indole: {reactants_with_indole}")

                # Check for N-sulfonylation (can occur at any depth >= 1)
                if not n_sulfonylation_detected and depth >= 1:
                    # Check for sulfonamide formation using sulfonyl halide
                    if checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ) or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    ):
                        print(f"Detected N-sulfonylation reaction at depth {depth}")
                        n_sulfonylation_detected = True
                    # Check for N-sulfonylation of indole or other heterocycles
                    elif checker.check_ring("indole", product_smiles) and any(
                        checker.check_fg("Sulfonyl halide", r) for r in reactants_smiles
                    ):
                        print(f"Detected N-sulfonylation of indole at depth {depth}")
                        n_sulfonylation_detected = True
                    # Check for general N-sulfonylation pattern
                    elif "S(=O)(=O)N" in product_smiles and any(
                        checker.check_fg("Sulfonyl halide", r) for r in reactants_smiles
                    ):
                        print(f"Detected general N-sulfonylation at depth {depth}")
                        n_sulfonylation_detected = True
                    else:
                        print(f"No N-sulfonylation detected at depth {depth}")
                        print(
                            f"Product has sulfonamide: {checker.check_fg('Sulfonamide', product_smiles)}"
                        )
                        print(f"Product has S(=O)(=O)N pattern: {'S(=O)(=O)N' in product_smiles}")
                        reactants_with_sulfonyl = [
                            r for r in reactants_smiles if checker.check_fg("Sulfonyl halide", r)
                        ]
                        print(f"Reactants with sulfonyl halide: {reactants_with_sulfonyl}")

                # Check for benzylic modification (can occur at any depth, but typically late stage)
                if not benzylic_modification_detected:
                    # Check for benzylic bromination
                    if "CH2Br" in product_smiles and not any(
                        "CH2Br" in r for r in reactants_smiles
                    ):
                        print(f"Detected benzylic bromination at depth {depth}")
                        benzylic_modification_detected = True
                    # Check for Wohl-Ziegler bromination
                    elif checker.check_reaction(
                        "Wohl-Ziegler bromination benzyl primary", rsmi
                    ) or checker.check_reaction("Wohl-Ziegler bromination benzyl secondary", rsmi):
                        print(f"Detected Wohl-Ziegler bromination at depth {depth}")
                        benzylic_modification_detected = True
                    # Check for other benzylic modifications
                    elif checker.check_fg("Primary halide", product_smiles) and not any(
                        checker.check_fg("Primary halide", r) for r in reactants_smiles
                    ):
                        print(f"Detected benzylic halogenation at depth {depth}")
                        benzylic_modification_detected = True
                    else:
                        print(f"No benzylic modification detected at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the complete strategy is detected
    strategy_detected = (
        indole_formation_detected and n_sulfonylation_detected and benzylic_modification_detected
    )

    print(f"Indole formation detected: {indole_formation_detected}")
    print(f"N-sulfonylation detected: {n_sulfonylation_detected}")
    print(f"Benzylic modification detected: {benzylic_modification_detected}")
    print(f"Complete strategy detected: {strategy_detected}")

    return strategy_detected
