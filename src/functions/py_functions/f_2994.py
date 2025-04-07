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
    Detects if the synthesis follows a linear strategy with sequential functional group transformations.
    """
    # Track synthesis metrics
    branch_count = 0
    max_depth = 0
    functional_group_changes = 0

    # List of functional groups to check
    functional_groups = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Carboxylic acid",
        "Aldehyde",
        "Ketone",
        "Alkyne",
        "Alkene",
        "Ester",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Nitrile",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Phenol",
        "Ether",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Boronic acid",
        "Boronic ester",
    ]

    # Common reaction types that indicate functional group transformations
    reaction_types = [
        "Oxidation of aldehydes to carboxylic acids",
        "Reduction of aldehydes and ketones to alcohols",
        "Esterification of Carboxylic Acids",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Williamson Ether Synthesis",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Reductive amination with aldehyde",
        "Reductive amination with ketone",
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of alcohol to aldehyde",
        "Reduction of carboxylic acid to primary alcohol",
        "Reduction of ester to primary alcohol",
        "Reduction of nitrile to amine",
    ]

    # Track the sequence of functional group transformations
    transformation_sequence = []

    def dfs_traverse(node, current_depth=0, path_id="0"):
        nonlocal branch_count, max_depth, functional_group_changes

        # Update max depth
        max_depth = max(max_depth, current_depth)

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check for functional group transformations
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Count complex reactants (non-reagents)
                    complex_reactants = [
                        r
                        for r in reactants
                        if Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).GetNumHeavyAtoms() > 3
                    ]

                    # Check for functional group changes
                    for reactant in complex_reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            fg_changes_in_reaction = 0
                            transformation = []

                            for fg in functional_groups:
                                reactant_has_fg = checker.check_fg(fg, reactant)
                                product_has_fg = checker.check_fg(fg, product)

                                if reactant_has_fg and not product_has_fg:
                                    transformation.append(f"{fg} removed")
                                    fg_changes_in_reaction += 1
                                elif not reactant_has_fg and product_has_fg:
                                    transformation.append(f"{fg} added")
                                    fg_changes_in_reaction += 1

                            # Check for specific reaction types
                            for rxn_type in reaction_types:
                                if checker.check_reaction(rxn_type, rsmi):
                                    transformation.append(f"Reaction: {rxn_type}")
                                    fg_changes_in_reaction += 1
                                    print(
                                        f"Detected reaction type: {rxn_type} at depth {current_depth}, path {path_id}"
                                    )

                            if fg_changes_in_reaction > 0:
                                functional_group_changes += 1
                                transformation_sequence.append(
                                    (current_depth, path_id, transformation)
                                )
                                print(
                                    f"Functional group change at depth {current_depth}, path {path_id}: {transformation}"
                                )
            except Exception as e:
                print(f"Error processing SMILES in linear_synthesis_strategy: {e}")

        # Count children for branching analysis
        children = node.get("children", [])
        if len(children) > 1:
            branch_count += 1
            print(
                f"Detected branch at depth {current_depth} with {len(children)} children"
            )

        # Continue traversing
        for i, child in enumerate(children):
            # Create a new path ID for each branch
            new_path_id = path_id if len(children) <= 1 else f"{path_id}-{i}"
            dfs_traverse(child, current_depth + 1, new_path_id)

    # Start traversal
    dfs_traverse(route)

    # Group transformations by path
    path_transformations = {}
    for depth, path_id, transformation in transformation_sequence:
        if path_id not in path_transformations:
            path_transformations[path_id] = []
        path_transformations[path_id].append((depth, transformation))

    # Sort each path's transformations by depth
    for path_id in path_transformations:
        path_transformations[path_id].sort(key=lambda x: x[0])

    # Find the longest path with sequential transformations
    longest_sequential_path = []
    for path_id, transformations in path_transformations.items():
        if len(transformations) > len(longest_sequential_path):
            # Check if transformations are sequential
            sequential = True
            if len(transformations) >= 2:
                for i in range(1, len(transformations)):
                    # Allow gaps of up to 5 in depth for sequential transformations
                    if transformations[i][0] - transformations[i - 1][0] > 5:
                        sequential = False
                        break
            else:
                sequential = False

            if sequential:
                longest_sequential_path = transformations

    print(f"Branch count: {branch_count}")
    print(f"Functional group changes: {functional_group_changes}")
    print(f"Max depth: {max_depth}")
    print(f"Longest sequential path: {longest_sequential_path}")

    # A linear synthesis has sequential functional group changes
    if len(longest_sequential_path) >= 2:
        print("Detected linear synthesis strategy")
        return True

    # Alternative criteria: few branches but many functional group changes
    if branch_count <= 3 and functional_group_changes >= 4 and max_depth >= 4:
        print("Detected linear synthesis strategy (alternative criteria)")
        return True

    # Check if we have a deep synthesis with multiple functional group changes
    if max_depth >= 8 and functional_group_changes >= 3:
        print("Detected linear synthesis strategy (deep synthesis)")
        return True

    return False
