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
    Detects if a complex polycyclic core structure is preserved throughout the synthesis.

    This strategy checks if a polycyclic core (multiple fused rings) is maintained throughout
    the main synthetic pathway, ignoring starting materials and reagents.
    """
    # Track molecules in the main synthetic pathway
    main_pathway_mols = []

    # Function to check if a molecule has a polycyclic core and identify which ones
    def identify_polycyclic_cores(mol_smiles):
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            print(f"Could not parse molecule: {mol_smiles}")
            return []

        # Check for common polycyclic ring systems
        polycyclic_rings = [
            "naphthalene",
            "anthracene",
            "purine",
            "carbazole",
            "acridine",
            "phenothiazine",
            "phenoxazine",
            "dibenzofuran",
            "dibenzothiophene",
            "xanthene",
            "thioxanthene",
        ]

        found_cores = []
        for ring in polycyclic_rings:
            if checker.check_ring(ring, mol_smiles):
                print(f"Found polycyclic core ({ring}) in: {mol_smiles}")
                found_cores.append(ring)

        # If no predefined cores found, check for fused ring systems using RingInfo
        if not found_cores:
            ri = mol.GetRingInfo()
            if ri.NumRings() >= 2:
                # Check if rings share atoms (indicating fusion)
                rings = ri.AtomRings()
                for i in range(len(rings)):
                    for j in range(i + 1, len(rings)):
                        # Check if rings share any atoms (indicating fusion)
                        if set(rings[i]).intersection(set(rings[j])):
                            print(f"Found generic fused ring system in: {mol_smiles}")
                            found_cores.append("fused_rings")
                            break
                    if "fused_rings" in found_cores:
                        break

        return found_cores

    # DFS traversal with depth tracking
    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            # Skip explicitly marked starting materials
            if node.get("in_stock", False):
                print(f"Skipping starting material: {node['smiles']}")
                return

            # Add to main pathway if not a starting material
            cores = identify_polycyclic_cores(node["smiles"])
            main_pathway_mols.append(
                {"smiles": node["smiles"], "depth": depth, "polycyclic_cores": cores}
            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If no main pathway molecules found, return False
    if not main_pathway_mols:
        print("No main pathway molecules found")
        return False

    # Sort by depth (target molecule has lowest depth)
    main_pathway_mols.sort(key=lambda x: x["depth"])
    target_mol = main_pathway_mols[0]

    # Check if the target molecule has a polycyclic core
    if not target_mol["polycyclic_cores"]:
        print("Target molecule does not have a polycyclic core")
        return False

    # Check if at least one of the target's polycyclic cores is preserved in intermediates
    preserved_cores = []
    for core in target_mol["polycyclic_cores"]:
        # Check if this core appears in any intermediate
        for mol in main_pathway_mols[1:]:
            if core in mol["polycyclic_cores"]:
                preserved_cores.append(core)
                break

    if not preserved_cores:
        print("No polycyclic cores are preserved in intermediates")
        return False

    print(
        f"Found {len(preserved_cores)} polycyclic cores preserved throughout synthesis: {preserved_cores}"
    )
    return True
