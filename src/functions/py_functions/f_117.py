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
    Detects if the synthetic route preserves heterocyclic structures throughout the synthesis.
    """
    # List of heterocyclic rings to check
    heterocycle_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "trioxane",
        "dioxepane",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "porphyrin",
        "indazole",
        "benzotriazole",
    ]

    # Track heterocycles through the synthesis
    preserved_heterocycles = {}  # Maps heterocycle type to preservation status

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            current_heterocycles = {}

            # Check for all heterocycle types
            for ring_type in heterocycle_types:
                if checker.check_ring(ring_type, mol_smiles):
                    # Get atom indices for this ring
                    ring_indices = checker.get_ring_atom_indices(ring_type, mol_smiles)
                    if ring_indices:
                        current_heterocycles[ring_type] = ring_indices
                        print(
                            f"Found {ring_type} in molecule at depth {depth}: {mol_smiles}"
                        )

            # If this is the final product (depth 0), initialize tracking
            if depth == 0:
                for ring_type, indices in current_heterocycles.items():
                    preserved_heterocycles[ring_type] = {
                        "in_final": True,
                        "in_starting": False,
                        "indices": indices,
                    }
                print(
                    f"Final product heterocycles: {list(current_heterocycles.keys())}"
                )

            # If this is a starting material, check if it contains heterocycles from final product
            elif node.get("in_stock", False):
                for ring_type, indices in current_heterocycles.items():
                    if ring_type in preserved_heterocycles:
                        preserved_heterocycles[ring_type]["in_starting"] = True
                        print(f"Starting material contains {ring_type}: {mol_smiles}")

        # Process children (reactants in retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path + [node])

    # Traverse the route
    dfs_traverse(route)

    # Check if any heterocycle is preserved from starting material to final product
    for ring_type, status in preserved_heterocycles.items():
        if status["in_final"] and status["in_starting"]:
            print(f"Heterocycle {ring_type} is preserved throughout the synthesis")
            return True

    print("No heterocycles preserved throughout the synthesis")
    return False
