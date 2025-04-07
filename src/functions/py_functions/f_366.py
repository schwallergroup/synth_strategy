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
    This function detects if the synthesis involves coupling between thiophene and pyridine rings.
    """
    has_thiophene_pyridine_coupling = False

    def dfs_traverse(node):
        nonlocal has_thiophene_pyridine_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if reactants contain thiophene and pyridine rings
                has_thiophene_reactant = any(
                    checker.check_ring("thiophene", r) for r in reactants_smiles
                )
                has_pyridine_reactant = any(
                    checker.check_ring("pyridine", r) for r in reactants_smiles
                )

                # Check if product has both thiophene and pyridine rings
                has_thiophene_product = checker.check_ring("thiophene", product_smiles)
                has_pyridine_product = checker.check_ring("pyridine", product_smiles)

                # Check if this is a coupling reaction
                is_coupling_reaction = (
                    checker.check_reaction("Suzuki", rsmi)
                    or checker.check_reaction("Stille", rsmi)
                    or checker.check_reaction("Negishi", rsmi)
                    or checker.check_reaction("Heck", rsmi)
                    or checker.check_reaction("Sonogashira", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                )

                # Check if thiophene and pyridine are in separate reactants
                separate_rings = False
                if has_thiophene_reactant and has_pyridine_reactant:
                    thiophene_reactants = [
                        r
                        for r in reactants_smiles
                        if checker.check_ring("thiophene", r)
                    ]
                    pyridine_reactants = [
                        r for r in reactants_smiles if checker.check_ring("pyridine", r)
                    ]

                    # Check if there's at least one reactant with thiophene but no pyridine
                    # and at least one with pyridine but no thiophene
                    thiophene_only = any(
                        not checker.check_ring("pyridine", r)
                        for r in thiophene_reactants
                    )
                    pyridine_only = any(
                        not checker.check_ring("thiophene", r)
                        for r in pyridine_reactants
                    )

                    separate_rings = thiophene_only and pyridine_only

                # Check if rings were already connected in any reactant
                rings_already_connected = False
                for reactant in reactants_smiles:
                    if checker.check_ring("thiophene", reactant) and checker.check_ring(
                        "pyridine", reactant
                    ):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            thiophene_indices = checker.get_ring_atom_indices(
                                "thiophene", reactant
                            )
                            pyridine_indices = checker.get_ring_atom_indices(
                                "pyridine", reactant
                            )

                            if (
                                thiophene_indices
                                and pyridine_indices
                                and len(thiophene_indices) > 0
                                and len(pyridine_indices) > 0
                            ):
                                if are_rings_connected(
                                    reactant_mol, thiophene_indices, pyridine_indices
                                ):
                                    rings_already_connected = True
                                    break

                # Verify that the product has both rings and they're connected
                if (
                    has_thiophene_product
                    and has_pyridine_product
                    and is_coupling_reaction
                    and separate_rings
                    and not rings_already_connected
                ):
                    # Create RDKit mol object for the product to check connectivity
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # Get atom indices for thiophene and pyridine rings in the product
                        thiophene_indices = checker.get_ring_atom_indices(
                            "thiophene", product_smiles
                        )
                        pyridine_indices = checker.get_ring_atom_indices(
                            "pyridine", product_smiles
                        )

                        if (
                            thiophene_indices
                            and pyridine_indices
                            and len(thiophene_indices) > 0
                            and len(pyridine_indices) > 0
                        ):
                            if are_rings_connected(
                                product_mol, thiophene_indices, pyridine_indices
                            ):
                                print(
                                    f"Found thiophene-pyridine coupling in reaction: {rsmi}"
                                )
                                has_thiophene_pyridine_coupling = True
            except Exception as e:
                print(f"Error processing reaction SMILES {rsmi}: {str(e)}")

        for child in node.get("children", []):
            dfs_traverse(child)

    def are_rings_connected(mol, ring1_indices, ring2_indices):
        """Check if two rings are connected in a molecule, either directly or through a linker."""
        # Flatten the nested tuples to get all atom indices for each ring
        ring1_atoms = set(
            [atom for match in ring1_indices for tup in match for atom in tup]
        )
        ring2_atoms = set(
            [atom for match in ring2_indices for tup in match for atom in tup]
        )

        # First check for direct bonds
        for atom1 in ring1_atoms:
            for atom2 in ring2_atoms:
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                if bond is not None:
                    return True

        # If no direct bonds, check for paths between rings
        # Use BFS to find paths between any atom in ring1 and any atom in ring2
        for start_atom in ring1_atoms:
            visited = set([start_atom])
            queue = deque([(start_atom, 0)])  # (atom_idx, distance)

            while queue:
                current_atom, distance = queue.popleft()

                # Limit path length to avoid considering very distant connections
                if distance > 3:  # Allow up to 3-atom linkers
                    continue

                # Get neighbors of current atom
                atom = mol.GetAtomWithIdx(current_atom)
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()

                    # Skip already visited atoms
                    if neighbor_idx in visited:
                        continue

                    # If neighbor is in ring2, we found a path
                    if neighbor_idx in ring2_atoms:
                        return True

                    # Otherwise, add to queue for further exploration
                    visited.add(neighbor_idx)
                    queue.append((neighbor_idx, distance + 1))

        return False

    dfs_traverse(route)
    return has_thiophene_pyridine_coupling
