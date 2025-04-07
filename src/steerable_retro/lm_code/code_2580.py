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
    This function detects if the synthesis uses a piperazine scaffold as a linker
    between two aromatic fragments.
    """
    piperazine_found = False
    final_product_has_piperazine = False
    connects_two_fragments = False

    def has_two_aromatic_fragments(mol_smiles, piperazine_indices=None):
        """Check if piperazine connects two aromatic fragments"""
        if piperazine_indices is None:
            piperazine_indices = checker.get_ring_atom_indices("piperazine", mol_smiles)
            if not piperazine_indices:
                print(f"No valid piperazine indices found: {piperazine_indices}")
                return False

        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            print(f"Could not create molecule from SMILES: {mol_smiles}")
            return False

        # Flatten the piperazine indices
        piperazine_atoms = set()
        try:
            # Handle different possible formats of piperazine_indices
            if isinstance(piperazine_indices, list) and len(piperazine_indices) > 0:
                if isinstance(piperazine_indices[0], int):
                    # Format is a simple list of integers
                    for atom_idx in piperazine_indices:
                        piperazine_atoms.add(atom_idx)
                elif isinstance(piperazine_indices[0], tuple):
                    # Format is a list of tuples
                    for ring_indices in piperazine_indices:
                        for atom_idx in ring_indices:
                            piperazine_atoms.add(atom_idx)
                else:
                    # Try to handle as a list of lists
                    for ring_indices in piperazine_indices:
                        if isinstance(ring_indices, (list, tuple)):
                            for atom_idx in ring_indices:
                                piperazine_atoms.add(atom_idx)
                        else:
                            piperazine_atoms.add(ring_indices)
            elif isinstance(piperazine_indices, tuple):
                # Format is a tuple
                for atom_idx in piperazine_indices:
                    piperazine_atoms.add(atom_idx)
        except (TypeError, IndexError) as e:
            print(f"Error processing piperazine indices: {piperazine_indices}, Error: {e}")
            # Instead of returning False, let's try a different approach
            # Let's get all atoms in piperazine rings directly from the molecule
            piperazine_atoms = set()
            patt = Chem.MolFromSmarts("N1CCN(C)CC1")  # Basic piperazine pattern
            if patt:
                matches = mol.GetSubstructMatches(patt)
                for match in matches:
                    for atom_idx in match:
                        piperazine_atoms.add(atom_idx)

            if not piperazine_atoms:
                return False

        if not piperazine_atoms:
            print("No piperazine atoms found after processing indices")
            return False

        # Find atoms connected to piperazine
        connected_atoms = set()
        for atom_idx in piperazine_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in piperazine_atoms:
                    connected_atoms.add(neighbor.GetIdx())

        # Check if connected atoms belong to aromatic fragments
        aromatic_fragments = 0
        visited = set()

        for start_atom in connected_atoms:
            if start_atom in visited:
                continue

            # BFS to find connected aromatic fragment
            queue = [start_atom]
            fragment = set()
            has_aromatic = False

            while queue:
                atom_idx = queue.pop(0)
                if atom_idx in visited or atom_idx in piperazine_atoms:
                    continue

                visited.add(atom_idx)
                fragment.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)

                if atom.GetIsAromatic():
                    has_aromatic = True

                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx not in visited and neighbor_idx not in piperazine_atoms:
                        queue.append(neighbor_idx)

            if has_aromatic and fragment:
                aromatic_fragments += 1

        print(f"Found {aromatic_fragments} aromatic fragments connected to piperazine")
        return aromatic_fragments >= 2

    def dfs_traverse(node):
        nonlocal piperazine_found, final_product_has_piperazine, connects_two_fragments

        if node["type"] == "mol":
            # Check if molecule contains piperazine
            if checker.check_ring("piperazine", node["smiles"]):
                piperazine_found = True
                print(f"Found piperazine in molecule: {node['smiles']}")

                # If this is the final product (root node)
                if node == route:
                    final_product_has_piperazine = True
                    print(f"Final product has piperazine: {node['smiles']}")

                    # Check if piperazine connects two aromatic fragments
                    if has_two_aromatic_fragments(node["smiles"]):
                        connects_two_fragments = True
                        print("Piperazine connects two aromatic fragments in final product")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reaction forms piperazine or connects fragments through piperazine
                product_has_piperazine = checker.check_ring("piperazine", product)

                if product_has_piperazine:
                    # Check if piperazine is newly formed or was already in reactants
                    reactants_with_piperazine = [
                        r for r in reactants if checker.check_ring("piperazine", r)
                    ]

                    if len(reactants) >= 2 and (len(reactants_with_piperazine) < len(reactants)):
                        # Either piperazine is newly formed or it's connecting to a new fragment
                        print(f"Reaction potentially forms or modifies piperazine: {rsmi}")

                        # Check if piperazine connects two aromatic fragments in the product
                        if has_two_aromatic_fragments(product):
                            connects_two_fragments = True
                            print(
                                f"Found reaction connecting aromatic fragments through piperazine: {rsmi}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    result = piperazine_found and final_product_has_piperazine and connects_two_fragments
    print(f"Piperazine linker strategy detected: {result}")
    return result
