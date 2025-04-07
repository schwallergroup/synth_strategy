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
    This function detects if the synthetic route involves a late-stage ring disconnection
    (early in retrosynthesis, which means low depth).
    """
    ring_disconnection_depth = None

    # List of common ring types to check
    ring_types = [
        "benzene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "piperidine",
        "tetrahydrofuran",
        "tetrahydropyran",
    ]

    def dfs_traverse(node, current_depth=0):
        nonlocal ring_disconnection_depth

        if node["type"] == "reaction":
            try:
                # Get depth from metadata or use current_depth
                depth = int(node["metadata"].get("depth", current_depth))
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for specific rings in product
                product_rings = {}
                for ring_type in ring_types:
                    if checker.check_ring(ring_type, product_smiles):
                        product_rings[ring_type] = True
                        print(f"Found {ring_type} in product")

                # Check for rings in reactants
                reactant_rings = {}
                for reactant in reactants_smiles:
                    for ring_type in ring_types:
                        if checker.check_ring(ring_type, reactant):
                            reactant_rings[ring_type] = reactant_rings.get(ring_type, 0) + 1
                            print(f"Found {ring_type} in reactant: {reactant}")

                # Check if any ring in product is not in reactants (ring disconnection)
                for ring_type in product_rings:
                    if ring_type not in reactant_rings:
                        print(f"Ring disconnection detected: {ring_type} at depth {depth}")
                        if ring_disconnection_depth is None or depth < ring_disconnection_depth:
                            ring_disconnection_depth = depth

                # Also check total ring count as a fallback
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactants_ring_count = 0
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        reactants_ring_count += mol.GetRingInfo().NumRings()

                if product_mol:
                    product_ring_count = product_mol.GetRingInfo().NumRings()
                    if product_ring_count > reactants_ring_count:
                        print(
                            f"Ring count disconnection at depth {depth}: Product has {product_ring_count} rings, reactants have {reactants_ring_count} rings"
                        )
                        if ring_disconnection_depth is None or depth < ring_disconnection_depth:
                            ring_disconnection_depth = depth

            except Exception as e:
                print(f"Error in ring disconnection detection: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it happens at depth 0, 1, or 2
    result = ring_disconnection_depth is not None and ring_disconnection_depth <= 2
    print(f"Late-stage ring disconnection detected: {result} (depth: {ring_disconnection_depth})")

    return result
