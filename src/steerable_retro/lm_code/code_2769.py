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
    Detects if an iodine substituent on a heterocycle is maintained throughout the synthesis.
    """
    # List of heterocyclic rings to check
    heterocycles = [
        "pyrazole",
        "pyrrole",
        "imidazole",
        "triazole",
        "tetrazole",
        "pyridine",
        "pyrimidine",
        "pyridazine",
        "pyrazine",
        "indole",
        "benzimidazole",
        "benzotriazole",
        "oxazole",
        "thiazole",
        "isoxazole",
    ]

    # Track molecules with iodine-substituted heterocycles
    iodine_heterocycle_molecules = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]

            # Check if molecule contains both a heterocycle and an iodine
            has_heterocycle = any(checker.check_ring(ring, mol_smiles) for ring in heterocycles)
            has_iodine = checker.check_fg("Aromatic halide", mol_smiles)

            if has_heterocycle and has_iodine:
                # Verify it's specifically an iodine (not just any halide)
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetAtomicNum() == 53:  # Iodine atomic number
                            # Check if iodine is attached to a heterocycle
                            for neighbor in atom.GetNeighbors():
                                if neighbor.IsInRing():
                                    iodine_heterocycle_molecules.append(
                                        {"smiles": mol_smiles, "depth": depth}
                                    )
                                    print(
                                        f"Found iodine-substituted heterocycle at depth {depth}: {mol_smiles}"
                                    )
                                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have iodine-heterocycles at different synthesis stages
    if len(iodine_heterocycle_molecules) < 2:
        return False

    # Sort by depth to check if maintained from early to late stages
    iodine_heterocycle_molecules.sort(key=lambda x: x["depth"])

    # Check if we have molecules at different depths (indicating maintenance through synthesis)
    depth_values = [mol["depth"] for mol in iodine_heterocycle_molecules]

    # If we have molecules at at least two different depths, and the difference between
    # min and max depth is at least 1, then the iodine-heterocycle is maintained
    if len(set(depth_values)) >= 2 and max(depth_values) - min(depth_values) >= 1:
        print(f"Iodine-heterocycle maintained through depths: {depth_values}")
        return True

    return False
