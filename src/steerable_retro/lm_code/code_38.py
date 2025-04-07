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
    Detects if the synthesis builds a complex molecule with multiple heterocyclic systems
    (e.g., pyrimidine, tetrazole, morpholine).
    """
    # List of heterocycles to check from the provided list
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

    # Get the target molecule (root node)
    target_mol_smiles = route["smiles"]

    # Set to track detected heterocycles
    detected_heterocycles = set()

    # Check for heterocycles in the target molecule
    for heterocycle in heterocycle_types:
        if checker.check_ring(heterocycle, target_mol_smiles):
            detected_heterocycles.add(heterocycle)
            print(f"Detected {heterocycle} heterocycle in target molecule")

    # If we don't find enough heterocycles in the target, check all molecules in the route
    if len(detected_heterocycles) < 3:

        def dfs_traverse(node):
            if node["type"] == "mol" and node.get("smiles"):
                mol_smiles = node["smiles"]

                # Skip if this is the target molecule (already checked)
                if mol_smiles == target_mol_smiles:
                    return

                # Check for each heterocycle type
                for heterocycle in heterocycle_types:
                    if heterocycle not in detected_heterocycles and checker.check_ring(
                        heterocycle, mol_smiles
                    ):
                        detected_heterocycles.add(heterocycle)
                        print(f"Detected {heterocycle} heterocycle in intermediate molecule")

            # Continue traversing
            for child in node.get("children", []):
                dfs_traverse(child)

        dfs_traverse(route)

    # Return True if we have 3 or more different heterocycle types
    print(f"Total heterocycles detected: {len(detected_heterocycles)}")
    print(f"Heterocycle types: {detected_heterocycles}")
    return len(detected_heterocycles) >= 3
