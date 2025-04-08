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
    This function detects if the synthesis involves formation of a heterocyclic ring system,
    specifically focusing on imidazopyridine-like structures.
    """
    ring_formation_detected = False

    # List of heterocyclic rings to check
    heterocyclic_rings = [
        "imidazole",
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "oxadiazole",
        "thiadiazole",
    ]

    def dfs_traverse(node):
        nonlocal ring_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check if any heterocyclic ring is formed in this reaction
                    reactant_rings_present = {}
                    product_rings_present = {}

                    # Check which heterocyclic rings are present in reactants
                    for ring_name in heterocyclic_rings:
                        reactant_rings_present[ring_name] = checker.check_ring(
                            ring_name, reactants_smiles
                        )

                    # Check which heterocyclic rings are present in product
                    for ring_name in heterocyclic_rings:
                        product_rings_present[ring_name] = checker.check_ring(
                            ring_name, product_smiles
                        )

                    # Check if any heterocyclic ring is in product but not in reactants
                    for ring_name in heterocyclic_rings:
                        if (
                            product_rings_present[ring_name]
                            and not reactant_rings_present[ring_name]
                        ):
                            print(f"Detected heterocyclic ring formation: {ring_name}")
                            ring_formation_detected = True
                            break

                    # Also check for specific ring formation reactions
                    ring_formation_reactions = [
                        "Paal-Knorr pyrrole synthesis",
                        "Formation of NOS Heterocycles",
                        "benzimidazole_derivatives_carboxylic-acid/ester",
                        "benzimidazole_derivatives_aldehyde",
                        "benzothiazole",
                        "benzoxazole_arom-aldehyde",
                        "benzoxazole_carboxylic-acid",
                        "thiazole",
                        "tetrazole_terminal",
                        "pyrazole",
                        "Fischer indole",
                        "indole",
                        "oxadiazole",
                    ]

                    for rxn_name in ring_formation_reactions:
                        if checker.check_reaction(rxn_name, rsmi):
                            print(f"Detected heterocyclic ring formation reaction: {rxn_name}")
                            ring_formation_detected = True
                            break

                except Exception as e:
                    print(f"Error in heterocyclic ring analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ring_formation_detected
