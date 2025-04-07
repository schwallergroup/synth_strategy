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
    Detects if the synthesis uses a convergent approach with late-stage ether formation
    between two complex fragments.
    """
    result = False

    def get_path_to_node(root, target_node, current_path=None):
        """Helper function to find the path from root to a specific node"""
        if current_path is None:
            current_path = []

        # Check if this is the target node (by reference)
        if root is target_node:
            return current_path

        # Check children
        for child in root.get("children", []):
            path = get_path_to_node(child, target_node, current_path + [root])
            if path is not None:
                return path

        return None

    def calculate_depth(node):
        """Calculate depth based on traversal from target molecule"""
        # For the root node (target molecule)
        if node == route:
            return 0

        # Try to find path from root to this node
        path = get_path_to_node(route, node)
        if path:
            return len(path)

        # Alternative: count steps from target molecule
        depth = 0
        current = route
        visited = set()

        def find_node_depth(current, target, current_depth=0):
            if current is target:
                return current_depth

            if id(current) in visited:
                return None

            visited.add(id(current))

            for child in current.get("children", []):
                result = find_node_depth(child, target, current_depth + 1)
                if result is not None:
                    return result

            return None

        node_depth = find_node_depth(route, node)
        return node_depth if node_depth is not None else 3  # Default to a reasonable depth

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Calculate depth based on position in synthesis tree
            depth = calculate_depth(node)
            print(f"Calculated depth: {depth}")

            # Check if this is a late-stage reaction (depth 0, 1, 2, or 3)
            if depth <= 3:
                print(f"Analyzing late-stage reaction at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # Check if we have multiple reactants (convergent)
                if len(reactants) >= 2:
                    print("Multiple reactants detected (potential convergent synthesis)")

                    # Check for ether formation using checker function
                    product_has_ether = checker.check_fg("Ether", product)
                    print(f"Product has ether: {product_has_ether}")

                    if product_has_ether:
                        # Check if this is a known ether formation reaction
                        is_ether_formation = (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                            or checker.check_reaction(
                                "Williamson Ether Synthesis (intra to epoxy)", rsmi
                            )
                            or checker.check_reaction("Alcohol to ether", rsmi)
                            or checker.check_reaction("{Williamson ether}", rsmi)
                        )
                        print(f"Known ether formation reaction: {is_ether_formation}")

                        if not is_ether_formation:
                            # If not a known ether formation reaction, check if ether is newly formed
                            ether_exists_in_reactants = any(
                                checker.check_fg("Ether", r) for r in reactants
                            )
                            print(f"Ether exists in reactants: {ether_exists_in_reactants}")

                            if not ether_exists_in_reactants:
                                # This is likely an ether formation reaction
                                is_ether_formation = True
                                print("Inferred ether formation based on functional group analysis")

                        # Additional check for ether formation by examining reactants and products
                        if not is_ether_formation and product_has_ether:
                            # Check for alcohol in reactants
                            alcohol_in_reactants = any(
                                checker.check_fg("Primary alcohol", r)
                                or checker.check_fg("Secondary alcohol", r)
                                or checker.check_fg("Tertiary alcohol", r)
                                or checker.check_fg("Aromatic alcohol", r)
                                or checker.check_fg("Phenol", r)
                                for r in reactants
                            )

                            # Check for leaving groups in reactants
                            leaving_group_in_reactants = any(
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                                or checker.check_fg("Tosylate", r)
                                or checker.check_fg("Mesylate", r)
                                or checker.check_fg("Triflate", r)
                                for r in reactants
                            )

                            if alcohol_in_reactants and leaving_group_in_reactants:
                                is_ether_formation = True
                                print(
                                    "Inferred ether formation based on alcohol and leaving group analysis"
                                )

                        if is_ether_formation:
                            print("Ether formation confirmed, checking reactant complexity")
                            # Check if reactants are complex (contain rings or have sufficient size)
                            complex_reactants = 0
                            for i, r in enumerate(reactants):
                                try:
                                    r_mol = Chem.MolFromSmiles(r)
                                    if r_mol:
                                        # Check if reactant has rings
                                        has_ring = False
                                        common_rings = [
                                            "benzene",
                                            "pyridine",
                                            "furan",
                                            "pyrrole",
                                            "thiophene",
                                            "cyclopentane",
                                            "cyclohexane",
                                            "pyran",
                                            "dioxane",
                                            "tetrahydrofuran",
                                            "tetrahydropyran",
                                            "piperidine",
                                            "morpholine",
                                            "indole",
                                            "quinoline",
                                            "naphthalene",
                                        ]

                                        for ring_name in common_rings:
                                            if checker.check_ring(ring_name, r):
                                                has_ring = True
                                                print(f"Reactant {i+1} contains {ring_name} ring")
                                                break

                                        # Check if reactant has sufficient size
                                        atom_count = r_mol.GetNumAtoms()
                                        print(f"Reactant {i+1} has {atom_count} atoms")

                                        if has_ring or atom_count >= 6:
                                            complex_reactants += 1
                                            print(f"Reactant {i+1} is complex")
                                except Exception as e:
                                    print(f"Error analyzing reactant {i+1}: {e}")

                            print(f"Complex reactants count: {complex_reactants}")
                            # If at least two complex reactants, this is a convergent synthesis
                            if complex_reactants >= 2:
                                print(
                                    f"SUCCESS: Detected convergent synthesis with late-stage ether formation at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                                result = True

        # Continue traversing
        for child in node.get("children", []):
            if not result:  # Stop traversing if we already found a match
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return result
