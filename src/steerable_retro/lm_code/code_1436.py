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
    This function detects a convergent synthesis with late-stage amide coupling
    between two complex fragments.
    """
    # Track if we found the pattern
    found_late_stage_amide = False

    # Track synthesis branches to identify convergent synthesis
    synthesis_branches = {}
    branch_counter = 0

    def assign_branches(node, current_branch=None):
        nonlocal branch_counter

        if node["type"] == "mol":
            # If this is a starting material, assign a new branch
            if node.get("in_stock", False):
                branch_counter += 1
                synthesis_branches[node["smiles"]] = branch_counter
                return branch_counter
            # Otherwise, inherit branch from children
            elif "children" in node and node["children"]:
                child_branches = []
                for child in node["children"]:
                    branch = assign_branches(child, current_branch)
                    if branch:
                        child_branches.append(branch)

                # If multiple branches converge, keep track of all of them
                if child_branches:
                    synthesis_branches[node["smiles"]] = child_branches
                    return child_branches[0]  # Return first branch for parent inheritance

            # Inherit branch from parent if not determined yet
            if current_branch and node["smiles"] not in synthesis_branches:
                synthesis_branches[node["smiles"]] = current_branch
                return current_branch

        elif node["type"] == "reaction":
            # For reaction nodes, process all children
            for child in node.get("children", []):
                assign_branches(child, current_branch)

        return current_branch

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_amide

        # Check reaction nodes at late stage (low depth)
        if (
            node["type"] == "reaction" and depth <= 3
        ):  # Increased depth limit to catch more potential late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction using the checker
                is_amide_coupling = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Schotten-Baumann_amide", rsmi
                    )  # Added additional reaction check
                )

                print(f"Is amide coupling reaction: {is_amide_coupling}")

                # If not detected as amide coupling, try to check manually
                if not is_amide_coupling:
                    # Split reactants
                    reactants = reactants_str.split(".")

                    # Check for acid/amine reactants and amide product
                    acid_found = any(
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Ester", r)
                        for r in reactants
                    )

                    amine_found = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )

                    amide_in_product = (
                        checker.check_fg("Primary amide", product_str)
                        or checker.check_fg("Secondary amide", product_str)
                        or checker.check_fg("Tertiary amide", product_str)
                    )

                    print(
                        f"Manual check - Acid: {acid_found}, Amine: {amine_found}, Amide in product: {amide_in_product}"
                    )

                    # If we have acid, amine, and amide, it's likely an amide coupling
                    if acid_found and amine_found and amide_in_product:
                        is_amide_coupling = True
                        print("Manually identified as amide coupling")

                if is_amide_coupling:
                    print(f"Found amide coupling reaction at depth {depth}")
                    # Split reactants
                    reactants = reactants_str.split(".")

                    # Need at least two reactants for convergent synthesis
                    if len(reactants) >= 2:
                        # Verify acid/acyl and amine in reactants, amide in product
                        acid_found = False
                        amine_found = False
                        acid_reactant = None
                        amine_reactant = None

                        # Check each reactant for acid/acyl or amine
                        for reactant in reactants:
                            try:
                                if (
                                    checker.check_fg("Carboxylic acid", reactant)
                                    or checker.check_fg("Acyl halide", reactant)
                                    or checker.check_fg("Ester", reactant)
                                ):
                                    acid_found = True
                                    acid_reactant = reactant
                                    print(f"Found acid/acyl reactant: {reactant}")
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Aniline", reactant)
                                ):
                                    amine_found = True
                                    amine_reactant = reactant
                                    print(f"Found amine reactant: {reactant}")
                            except Exception as e:
                                print(f"Error checking functional groups: {e}")
                                continue

                        # Check product for amide
                        try:
                            amide_in_product = (
                                checker.check_fg("Primary amide", product_str)
                                or checker.check_fg("Secondary amide", product_str)
                                or checker.check_fg("Tertiary amide", product_str)
                            )
                            print(f"Amide in product: {amide_in_product}")
                        except Exception as e:
                            print(f"Error checking product for amide: {e}")
                            amide_in_product = False

                        # If we have acid/acyl, amine, and amide formation, confirm it's an amide coupling
                        if acid_found and amine_found and amide_in_product:
                            print("Confirmed acid/acyl and amine reactants with amide product")

                            # Check if both reactants are complex
                            complex_reactants = []
                            for reactant in [acid_reactant, amine_reactant]:
                                if not reactant:
                                    continue
                                try:
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        # Define complexity as having rings or being large
                                        num_rings = mol.GetRingInfo().NumRings()
                                        num_atoms = mol.GetNumAtoms()
                                        num_hetero = sum(
                                            1
                                            for atom in mol.GetAtoms()
                                            if atom.GetAtomicNum() not in [1, 6]
                                        )

                                        print(
                                            f"Reactant complexity: {reactant} (rings: {num_rings}, atoms: {num_atoms}, heteroatoms: {num_hetero})"
                                        )

                                        # More flexible complexity check
                                        if (
                                            (num_rings >= 1 and num_atoms >= 7)
                                            or (num_atoms >= 10)
                                            or (num_rings >= 1 and num_hetero >= 1)
                                            or (num_hetero >= 3)
                                        ):
                                            complex_reactants.append(reactant)
                                            print(f"Complex reactant found: {reactant}")
                                except Exception as e:
                                    print(f"Error checking reactant complexity: {e}")
                                    continue

                            # Check if we have at least 2 complex reactants
                            if len(complex_reactants) >= 2:
                                print("Found at least 2 complex reactants")

                                # Check if reactants come from different branches (convergent synthesis)
                                if (
                                    acid_reactant in synthesis_branches
                                    and amine_reactant in synthesis_branches
                                ):

                                    acid_branch = synthesis_branches[acid_reactant]
                                    amine_branch = synthesis_branches[amine_reactant]

                                    print(
                                        f"Acid branch: {acid_branch}, Amine branch: {amine_branch}"
                                    )

                                    # Check if branches are different
                                    if isinstance(acid_branch, list) or isinstance(
                                        amine_branch, list
                                    ):
                                        # If either is a list, check for non-overlapping branches
                                        acid_branches = (
                                            acid_branch
                                            if isinstance(acid_branch, list)
                                            else [acid_branch]
                                        )
                                        amine_branches = (
                                            amine_branch
                                            if isinstance(amine_branch, list)
                                            else [amine_branch]
                                        )

                                        # Check if there's at least one branch that's different
                                        if set(acid_branches) != set(amine_branches) and not (
                                            set(acid_branches).issubset(set(amine_branches))
                                            or set(amine_branches).issubset(set(acid_branches))
                                        ):
                                            print(
                                                f"Found convergent synthesis: acid branch {acid_branches}, amine branch {amine_branches}"
                                            )
                                            found_late_stage_amide = True
                                    elif acid_branch != amine_branch:
                                        print(
                                            f"Found convergent synthesis: acid branch {acid_branch}, amine branch {amine_branch}"
                                        )
                                        found_late_stage_amide = True
                                else:
                                    print(
                                        f"Reactants not found in synthesis branches: acid_reactant in branches: {acid_reactant in synthesis_branches}, amine_reactant in branches: {amine_reactant in synthesis_branches}"
                                    )
                            else:
                                print(
                                    f"Not enough complex reactants: found {len(complex_reactants)}"
                                )
                        else:
                            print(
                                f"Missing required functional groups: acid_found={acid_found}, amine_found={amine_found}, amide_in_product={amide_in_product}"
                            )
                    else:
                        print(f"Not enough reactants: found {len(reactants)}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First assign branches to track convergent synthesis
    assign_branches(route)
    print(f"Synthesis branches: {synthesis_branches}")

    # Then start traversal to find late-stage amide coupling
    dfs_traverse(route)

    # If we didn't find a late-stage amide coupling, try a more lenient approach
    if not found_late_stage_amide:
        print("\nTrying more lenient approach...")
        # Check if we have the right starting materials for a potential convergent synthesis
        starting_materials = [
            smiles for smiles, branch in synthesis_branches.items() if isinstance(branch, int)
        ]

        # Check if we have at least one carboxylic acid and one amine in starting materials
        acid_starting = any(checker.check_fg("Carboxylic acid", sm) for sm in starting_materials)
        amine_starting = any(
            checker.check_fg("Primary amine", sm)
            or checker.check_fg("Secondary amine", sm)
            or checker.check_fg("Aniline", sm)
            for sm in starting_materials
        )

        print(f"Starting materials contain acid: {acid_starting}, amine: {amine_starting}")

        # If we have both acid and amine starting materials from different branches, and the target molecule has an amide
        if acid_starting and amine_starting and len(set(synthesis_branches.values())) >= 2:
            # Check if the final product has an amide
            if route["type"] == "mol" and (
                checker.check_fg("Primary amide", route["smiles"])
                or checker.check_fg("Secondary amide", route["smiles"])
                or checker.check_fg("Tertiary amide", route["smiles"])
            ):
                print(
                    "Found amide in final product with acid and amine starting materials from different branches"
                )
                found_late_stage_amide = True

    return found_late_stage_amide
