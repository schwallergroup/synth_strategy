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
    This function detects if the synthesis follows a linear strategy (vs convergent).

    A linear synthesis is characterized by:
    1. Each reaction step has at most one complex reactant
    2. The synthesis tree has a predominantly linear structure

    Complex reactants are defined as:
    - Having more than 10 heavy atoms
    - Not being commercially available (in_stock)
    - Requiring multiple synthesis steps
    """
    is_linear = True
    convergent_steps_found = 0

    def assess_complexity(mol_smiles):
        """Assess the complexity of a molecule based on multiple factors"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return 0

        # Base complexity on heavy atom count
        complexity = mol.GetNumHeavyAtoms()

        # Add complexity for rings
        ring_info = mol.GetRingInfo()
        complexity += ring_info.NumRings() * 2

        # Add complexity for certain functional groups
        for fg in ["Ester", "Amide", "Carboxylic acid", "Nitro group", "Nitrile"]:
            if checker.check_fg(fg, mol_smiles):
                complexity += 2

        # Add complexity for complex ring systems
        for ring in ["indole", "quinoline", "naphthalene", "benzimidazole"]:
            if checker.check_ring(ring, mol_smiles):
                complexity += 5

        return complexity

    def is_complex_reactant(mol_node, complexity_threshold=15):
        """Determine if a molecule node represents a complex reactant"""
        # Check if it's in stock (commercially available)
        if mol_node.get("in_stock", False):
            return False

        # Check complexity based on structure
        complexity = assess_complexity(mol_node["smiles"])

        # Check if it requires further synthesis
        requires_synthesis = len(mol_node.get("children", [])) > 0

        return complexity > complexity_threshold and requires_synthesis

    def dfs_traverse(node, depth=0, path=None):
        nonlocal is_linear, convergent_steps_found

        if path is None:
            path = []

        # Add current node to path
        current_path = path + [node]

        # Check reaction nodes for convergent steps
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Get the reactant molecules from children
            reactant_nodes = [child for child in node.get("children", []) if child["type"] == "mol"]

            # Count complex reactants
            complex_reactant_count = sum(
                1 for r_node in reactant_nodes if is_complex_reactant(r_node)
            )

            # If more than one complex reactant, it's likely convergent
            if complex_reactant_count > 1:
                # Late-stage convergent steps are stronger indicators
                if depth < 3:  # Adjust threshold as needed
                    print(
                        f"Late-stage convergent step detected with {complex_reactant_count} complex reactants at depth {depth}"
                    )
                    is_linear = False
                    convergent_steps_found += 1
                else:
                    # Early-stage convergent steps might still be part of a linear strategy
                    print(
                        f"Early-stage convergent step detected with {complex_reactant_count} complex reactants at depth {depth}"
                    )
                    convergent_steps_found += 1
                    # Only mark as non-linear if we find multiple such steps
                    if convergent_steps_found > 2:
                        is_linear = False

        # Check for branching patterns typical in convergent synthesis
        if node["type"] == "mol" and not node.get("in_stock", False):
            # Analyze branches from this molecule
            reaction_branches = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]

            # Count branches that lead to complex molecules
            complex_branches = 0
            for rxn_branch in reaction_branches:
                complex_reactants_in_branch = 0
                for mol_child in [c for c in rxn_branch.get("children", []) if c["type"] == "mol"]:
                    if is_complex_reactant(mol_child):
                        complex_reactants_in_branch += 1

                if complex_reactants_in_branch > 0:
                    complex_branches += 1

            # Multiple complex branches indicate convergent synthesis
            if complex_branches > 1:
                print(
                    f"Convergent structure detected with {complex_branches} complex branches at depth {depth}"
                )
                is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    dfs_traverse(route)

    # Final check: if we found convergent steps but not enough to mark as non-linear
    if convergent_steps_found > 0 and is_linear:
        print(
            f"Found {convergent_steps_found} convergent steps, but still considering linear overall"
        )

    return is_linear
