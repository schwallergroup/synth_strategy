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
    Detects late-stage convergent amide coupling (in the first half of the synthesis).
    """
    found_late_amide_coupling = False
    max_depth = 0

    # First pass to determine the maximum depth of the route
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)
        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)
    first_half_threshold = max_depth // 2 if max_depth > 0 else 0

    print(f"Maximum depth: {max_depth}, First half threshold: {first_half_threshold}")

    def dfs_traverse(node, depth=0):
        nonlocal found_late_amide_coupling

        if (
            node["type"] == "reaction" and depth <= first_half_threshold
        ):  # Only consider reactions in first half
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                is_amide_coupling = False
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_coupling = True
                        print(f"Found amide coupling reaction: {reaction_type}")
                        break

                # If not found through reaction types, check for functional group changes
                if not is_amide_coupling:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for amide formation
                    has_amide_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    # Check for acid/acid derivative and amine in reactants
                    has_acid_derivative = False
                    has_amine = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_acid_derivative = True
                            print(f"Found acid derivative in reactant: {reactant}")

                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine = True
                            print(f"Found amine in reactant: {reactant}")

                    if has_acid_derivative and has_amine and has_amide_product:
                        is_amide_coupling = True
                        print("Detected amide coupling through functional group analysis")

                if is_amide_coupling:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if this is a convergent step (multiple complex fragments)
                    if len(reactants) >= 2:
                        # Check complexity of fragments
                        complex_fragments = 0
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Consider a fragment complex if it has >7 atoms and contains at least one ring
                                # or has >10 atoms with multiple functional groups
                                atom_count = reactant_mol.GetNumAtoms()
                                has_ring = reactant_mol.GetRingInfo().NumRings() > 0

                                # Count functional groups
                                fg_count = 0
                                for fg in [
                                    "Carboxylic acid",
                                    "Ester",
                                    "Alcohol",
                                    "Amine",
                                    "Amide",
                                    "Nitrile",
                                    "Halide",
                                ]:
                                    if checker.check_fg(fg, reactant):
                                        fg_count += 1

                                is_complex = (atom_count > 7 and has_ring) or (
                                    atom_count > 10 and fg_count >= 2
                                )

                                if is_complex:
                                    complex_fragments += 1
                                    print(
                                        f"Found complex fragment: {reactant} (atoms: {atom_count}, has ring: {has_ring}, fg_count: {fg_count})"
                                    )

                        if complex_fragments >= 2:
                            found_late_amide_coupling = True
                            print("Found late-stage convergent amide coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_amide_coupling
