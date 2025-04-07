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
    This function detects if the route contains C-C bond formation between aromatic systems
    via a tertiary alcohol intermediate.
    """
    found_cc_coupling = False

    def dfs_traverse(node):
        nonlocal found_cc_coupling

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for coupling between aromatic halide and ketone
            if len(reactants) >= 2:
                # Check for aromatic halide in reactants
                has_aromatic_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                )
                print(f"Has aromatic halide: {has_aromatic_halide}")

                # Check for ketone in reactants
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                print(f"Has ketone: {has_ketone}")

                # Check for tertiary alcohol in product
                has_tertiary_alcohol = checker.check_fg("Tertiary alcohol", product)
                print(f"Has tertiary alcohol: {has_tertiary_alcohol}")

                # Check if the reaction is a Grignard or similar coupling
                is_grignard_reaction = checker.check_reaction(
                    "Grignard from ketone to alcohol", rsmi
                )
                print(f"Is Grignard reaction: {is_grignard_reaction}")

                # Check for other relevant organometallic couplings
                organometallic_reactions = [
                    "Grignard from ketone to alcohol",
                    "Reaction of alkyl halides with organometallic coumpounds",
                    "Preparation of organolithium compounds",
                    "Aryllithium cross-coupling",
                    "Grignard with CO2 to carboxylic acid",
                    "Grignard from aldehyde to alcohol",
                    "Olefination of ketones with Grignard reagents",
                    "Olefination of aldehydes with Grignard reagents",
                ]

                has_organometallic_coupling = any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in organometallic_reactions
                )

                # Additional check for organometallic reagents in the reaction
                has_organometallic_reagent = (
                    "[Li]" in rsmi or "[Mg]" in rsmi or "Li" in rsmi or "Mg" in rsmi
                )
                print(f"Has organometallic coupling: {has_organometallic_coupling}")
                print(f"Has organometallic reagent: {has_organometallic_reagent}")

                # Check if product has at least two aromatic rings
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    aromatic_rings = 0
                    for ring in product_mol.GetRingInfo().AtomRings():
                        if all(
                            product_mol.GetAtomWithIdx(idx).GetIsAromatic()
                            for idx in ring
                        ):
                            aromatic_rings += 1
                    has_two_aromatic_rings = aromatic_rings >= 2
                    print(f"Number of aromatic rings in product: {aromatic_rings}")
                else:
                    has_two_aromatic_rings = False
                    print("Failed to parse product molecule")

                # Check if the reaction involves C-C bond formation
                if (
                    has_aromatic_halide
                    and has_tertiary_alcohol
                    and has_two_aromatic_rings
                    and (
                        is_grignard_reaction
                        or has_organometallic_coupling
                        or has_organometallic_reagent
                    )
                ):
                    # Additional check to ensure the tertiary alcohol carbon is connected to two aromatic systems
                    if product_mol:
                        # Find tertiary alcohol carbon
                        tertiary_alcohol_indices = checker.get_fg_atom_indices(
                            "Tertiary alcohol", product
                        )
                        if tertiary_alcohol_indices:
                            for indices_tuple in tertiary_alcohol_indices:
                                for indices in indices_tuple:
                                    # Find the carbon connected to OH in tertiary alcohol
                                    # In tertiary alcohol, find the carbon connected to oxygen
                                    carbon_idx = None
                                    oxygen_idx = None

                                    for idx in indices:
                                        atom = product_mol.GetAtomWithIdx(idx)
                                        if atom.GetSymbol() == "O":
                                            oxygen_idx = idx
                                            break

                                    if oxygen_idx is not None:
                                        oxygen_atom = product_mol.GetAtomWithIdx(
                                            oxygen_idx
                                        )
                                        for neighbor in oxygen_atom.GetNeighbors():
                                            if neighbor.GetSymbol() == "C":
                                                carbon_idx = neighbor.GetIdx()
                                                break

                                    if carbon_idx is None:
                                        continue

                                    # BFS to find aromatic systems within 3 bonds
                                    visited = set([carbon_idx])
                                    queue = deque(
                                        [(carbon_idx, 0)]
                                    )  # (atom_idx, distance)
                                    aromatic_systems = set()

                                    while queue:
                                        current_idx, distance = queue.popleft()
                                        current_atom = product_mol.GetAtomWithIdx(
                                            current_idx
                                        )

                                        # If we're at an aromatic atom, mark its ring system
                                        if (
                                            current_atom.GetIsAromatic()
                                            and distance > 0
                                        ):
                                            # Find which ring this atom belongs to
                                            for ring_idx, ring in enumerate(
                                                product_mol.GetRingInfo().AtomRings()
                                            ):
                                                if current_idx in ring and all(
                                                    product_mol.GetAtomWithIdx(
                                                        idx
                                                    ).GetIsAromatic()
                                                    for idx in ring
                                                ):
                                                    aromatic_systems.add(ring_idx)
                                                    break

                                        # Stop exploring this path if we're too far
                                        if distance >= 3:
                                            continue

                                        # Add neighbors to queue
                                        for neighbor in current_atom.GetNeighbors():
                                            neighbor_idx = neighbor.GetIdx()
                                            if neighbor_idx not in visited:
                                                visited.add(neighbor_idx)
                                                queue.append(
                                                    (neighbor_idx, distance + 1)
                                                )

                                    if len(aromatic_systems) >= 2:
                                        found_cc_coupling = True
                                        print(
                                            "Found C-C coupling between aromatic systems via tertiary alcohol"
                                        )
                                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_cc_coupling
