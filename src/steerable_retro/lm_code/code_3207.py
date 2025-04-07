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
    This function detects a convergent synthesis with late-stage amide coupling.
    It identifies if the final step combines two fragments via amide bond formation.
    """
    # Track if we found a convergent synthesis with amide coupling
    found_convergent_amide = False
    # Track branches and their reactants
    branch_reactants = {}
    # Track amide coupling reactions and their depths
    amide_couplings = []

    def dfs_traverse(node, depth=0, branch_id=0):
        nonlocal found_convergent_amide, branch_reactants, amide_couplings

        print(f"Examining node at depth {depth}, type: {node['type']}")

        # Check for reaction nodes at early stages (depth <= 1)
        if node["type"] == "reaction" and depth <= 1:
            print(f"Examining potential late-stage reaction at depth {depth}")

            # Check if this reaction has metadata with reaction SMILES
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reaction at depth {depth} has {len(reactants)} reactants")
                print(f"Reaction SMILES: {rsmi}")

                # For convergent synthesis, we need at least 2 reactants and 2 child nodes
                if len(reactants) >= 2 and len(node.get("children", [])) >= 2:
                    # Check for amide formation using checker functions
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )

                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    # Check if this is an amide coupling reaction using various reaction types
                    is_amide_coupling = (
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                        )
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                        )
                        or checker.check_reaction(
                            "Acyl chloride with secondary amine to amide", rsmi
                        )
                        or checker.check_reaction("Ester with primary amine to amide", rsmi)
                        or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    )

                    print(
                        f"Acid: {has_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}, Amine: {has_amine}, Amide: {has_amide}"
                    )
                    print(f"Is amide coupling reaction: {is_amide_coupling}")

                    # Consider it an amide coupling if either:
                    # 1. It matches a known amide coupling reaction type, or
                    # 2. It has the right functional groups and forms an amide
                    if is_amide_coupling or (
                        (has_acid or has_acyl_halide or has_ester) and has_amine and has_amide
                    ):
                        found_convergent_amide = True
                        amide_couplings.append((depth, rsmi))
                        print(f"Found amide coupling at depth {depth}!")

                        # Check if the reactants come from different branches
                        mol_nodes = [
                            child for child in node.get("children", []) if child["type"] == "mol"
                        ]
                        if len(mol_nodes) >= 2:
                            print(
                                f"Found at least 2 reactant molecules from potentially different branches"
                            )

        # Process children with appropriate branch IDs
        for i, child in enumerate(node.get("children", [])):
            # If this is a reaction node with multiple children, assign new branch IDs
            if node["type"] == "reaction" and len(node.get("children", [])) > 1:
                new_branch_id = branch_id + i
                print(f"Found branching at depth {depth}, assigning branch ID {new_branch_id}")
            else:
                new_branch_id = branch_id

            # Store molecule nodes by branch for later analysis
            if child["type"] == "mol" and "smiles" in child:
                if new_branch_id not in branch_reactants:
                    branch_reactants[new_branch_id] = []
                branch_reactants[new_branch_id].append(child["smiles"])

            dfs_traverse(child, depth + 1, new_branch_id)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final result: found_convergent_amide={found_convergent_amide}, branches={len(branch_reactants)}"
    )
    print(f"Amide couplings found: {amide_couplings}")

    # For a true convergent synthesis with late-stage amide coupling:
    # 1. We must have found an amide coupling reaction
    # 2. We must have at least 2 different branches (branch_reactants should have at least 2 keys)
    return found_convergent_amide and len(branch_reactants) >= 2
