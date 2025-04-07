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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    # Track the number of reactants at each step and branching points
    multi_reactant_steps = 0
    branching_points = 0

    def dfs_traverse(node):
        nonlocal multi_reactant_steps, branching_points

        # Check for branching in the synthesis tree
        if node["type"] == "mol":
            reaction_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "reaction"
            ]
            if len(reaction_children) > 1:
                print(f"Found branching point at molecule: {node['smiles']}")
                branching_points += 1

        # Check for convergent steps (multiple complex reactants)
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If a step has more than one reactant (excluding solvents, catalysts)
                if len(reactants) > 1:
                    # Count reactants with more than 10 heavy atoms as complex
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumHeavyAtoms() > 10:
                            # Don't count common activating agents as complex reactants
                            if not any(
                                checker.check_fg(fg, reactant)
                                for fg in [
                                    "Acyl halide",
                                    "Sulfonyl halide",
                                    "Triflate",
                                    "Mesylate",
                                    "Tosylate",
                                ]
                            ):
                                complex_reactants += 1

                    # Check if this is a convergent step, excluding common coupling reactions
                    if complex_reactants > 1:
                        # Some coupling reactions are commonly used in linear synthesis
                        common_coupling_rxns = [
                            "Suzuki coupling with boronic acids",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Sonogashira alkyne_aryl halide",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "Stille reaction_aryl",
                        ]

                        if not any(
                            checker.check_reaction(rxn, rsmi)
                            for rxn in common_coupling_rxns
                        ):
                            print(
                                f"Found convergent step with {complex_reactants} complex reactants: {rsmi}"
                            )
                            multi_reactant_steps += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If there are no convergent steps or branching points, it's a linear synthesis
    is_linear = multi_reactant_steps == 0 and branching_points == 0
    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print("Detected non-linear synthesis strategy")

    return is_linear
