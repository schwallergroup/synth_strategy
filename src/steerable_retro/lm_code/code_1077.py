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
    This function detects a synthesis pathway where a cross-coupling reaction
    is followed by a nucleophilic aromatic substitution on a heterocycle.
    """
    # Track reactions and their depths
    reactions_by_depth = {}

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "reaction":
            # Try to get depth from metadata, or use current_depth from DFS
            depth = node.get("metadata", {}).get("depth", current_depth)

            # Store reaction by depth
            if depth not in reactions_by_depth:
                reactions_by_depth[depth] = []

            reactions_by_depth[depth].append(node)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Found reactions at depths: {sorted(reactions_by_depth.keys())}")

    # Check for the pattern
    found_cross_coupling = False
    found_nucleophilic_substitution = False
    correct_sequence = False

    # Sort depths to ensure correct sequence checking
    depths = sorted(reactions_by_depth.keys())

    cross_coupling_depth = -1
    nucleophilic_subst_depth = -1

    for depth in depths:
        for reaction in reactions_by_depth[depth]:
            try:
                rsmi = reaction["metadata"]["rsmi"]
                print(f"Processing reaction at depth {depth}: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for cross-coupling using checker function
                is_cross_coupling = False

                # Check for Suzuki coupling
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                ):
                    print(f"Found Suzuki coupling at depth {depth}")
                    is_cross_coupling = True

                # Check for other cross-coupling reactions
                elif (
                    checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Heck terminal_vinyl", rsmi)
                    or checker.check_reaction("Heck_non-terminal_vinyl", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                ):
                    print(f"Found other cross-coupling reaction at depth {depth}")
                    is_cross_coupling = True

                # Manual check for cross-coupling pattern if checker fails
                if not is_cross_coupling:
                    # Check for boronic acid/ester and halide in reactants
                    has_boronic = any("B(O)" in r or "BO" in r for r in reactants)
                    has_halide = any(
                        ("Br" in r or "Cl" in r or "I" in r) and not ("N" in r and "Cl" in r)
                        for r in reactants
                    )

                    if has_boronic and has_halide:
                        print(f"Manually identified potential Suzuki coupling at depth {depth}")
                        is_cross_coupling = True

                if is_cross_coupling:
                    found_cross_coupling = True
                    cross_coupling_depth = depth

                # Check for nucleophilic aromatic substitution on a heterocycle
                is_nucleophilic_sub = False

                if (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                ):
                    print(f"Found nucleophilic substitution reaction at depth {depth}")
                    is_nucleophilic_sub = True

                # Manual check for nucleophilic substitution pattern
                if not is_nucleophilic_sub:
                    # Check for amine and halogenated heterocycle pattern
                    has_amine = any("NH2" in r or "NH" in r for r in reactants)
                    has_halogen_heterocycle = False

                    for r in reactants:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            # Check if molecule contains both a halogen and a heterocycle
                            has_halogen = "Cl" in r or "Br" in r or "I" in r or "F" in r
                            has_heterocycle = any(
                                atom.IsInRing() and atom.GetAtomicNum() in [7, 8, 16]
                                for atom in mol.GetAtoms()
                            )

                            if has_halogen and has_heterocycle:
                                has_halogen_heterocycle = True
                                break

                    if has_amine and has_halogen_heterocycle:
                        print(
                            f"Manually identified potential nucleophilic substitution at depth {depth}"
                        )
                        is_nucleophilic_sub = True

                if is_nucleophilic_sub:
                    # Check if the product contains a heterocycle
                    heterocycle_present = False
                    for ring_name in [
                        "pyridine",
                        "pyrimidine",
                        "pyrazine",
                        "pyridazine",
                        "triazole",
                        "tetrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "furan",
                        "thiophene",
                        "pyrrole",
                        "isoxazole",
                        "isothiazole",
                        "oxadiazole",
                        "thiadiazole",
                        "benzoxazole",
                        "benzothiazole",
                        "benzimidazole",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "purine",
                    ]:
                        if checker.check_ring(ring_name, product):
                            heterocycle_present = True
                            print(f"Found heterocycle {ring_name} in product")
                            break

                    # If no specific heterocycle was found, check for general heterocycle pattern
                    if not heterocycle_present:
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            heterocycle_present = any(
                                atom.IsInRing() and atom.GetAtomicNum() in [7, 8, 16]
                                for atom in mol.GetAtoms()
                            )
                            if heterocycle_present:
                                print("Found general heterocycle pattern in product")

                    if heterocycle_present:
                        print(
                            f"Found nucleophilic aromatic substitution on heterocycle at depth {depth}"
                        )
                        found_nucleophilic_substitution = True
                        nucleophilic_subst_depth = depth
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

    print(f"Cross-coupling found: {found_cross_coupling} at depth {cross_coupling_depth}")
    print(
        f"Nucleophilic substitution found: {found_nucleophilic_substitution} at depth {nucleophilic_subst_depth}"
    )

    # Check if the sequence is correct (cross-coupling followed by nucleophilic substitution)
    if found_cross_coupling and found_nucleophilic_substitution:
        # In retrosynthetic direction, nucleophilic substitution should be at a lower depth
        # than cross-coupling (meaning it happens later in the forward synthesis)
        if nucleophilic_subst_depth < cross_coupling_depth:
            print("Correct sequence: cross-coupling followed by nucleophilic substitution")
            correct_sequence = True
        else:
            print("Incorrect sequence: nucleophilic substitution does not follow cross-coupling")

    return correct_sequence
