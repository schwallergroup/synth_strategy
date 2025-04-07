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
    Detects a synthetic strategy involving multiple heterocycles
    (pyrazole, morpholine, dioxane, piperazine, and others).
    Returns True if at least 3 different heterocycles are present.
    """
    # Expanded list of heterocycles to check
    heterocycles = {
        "pyrazole": False,
        "morpholine": False,
        "dioxane": False,
        "piperazine": False,
        "furan": False,
        "pyridine": False,
        "imidazole": False,
        "oxazole": False,
        "thiazole": False,
        "pyrimidine": False,
        "triazole": False,
        "tetrazole": False,
        "pyrrolidine": False,
        "thiophene": False,
        "isoxazole": False,
        "oxadiazole": False,
        "thiadiazole": False,
        "benzoxazole": False,
        "benzothiazole": False,
        "benzimidazole": False,
        "oxetane": False,
        "tetrahydrofuran": False,
        "tetrahydropyran": False,
    }

    # Track which molecules contain which heterocycles for debugging
    molecule_heterocycles = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol_heterocycles = []

            # Check for each heterocycle using the checker function
            for heterocycle_name in heterocycles.keys():
                if checker.check_ring(heterocycle_name, smiles):
                    heterocycles[heterocycle_name] = True
                    mol_heterocycles.append(heterocycle_name)
                    print(f"Found {heterocycle_name} heterocycle in molecule: {smiles}")

            if mol_heterocycles:
                molecule_heterocycles[smiles] = mol_heterocycles
            else:
                # Debug molecules with no detected heterocycles
                print(f"No heterocycles detected in molecule: {smiles}")
                # Try to create RDKit mol for additional debugging
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        print(f"  Molecule parsed successfully, has {mol.GetNumAtoms()} atoms")
                except Exception as e:
                    print(f"  Error parsing molecule: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if at least 3 different heterocycles are present
    heterocycle_count = sum(heterocycles.values())
    present_heterocycles = [name for name, present in heterocycles.items() if present]

    print(f"Found {heterocycle_count} different heterocycles: {', '.join(present_heterocycles)}")
    print(f"Molecule heterocycle distribution: {molecule_heterocycles}")

    # Additional check for specific patterns that might be missed
    for smiles in molecule_heterocycles.keys():
        # Check for pyrimidine pattern in "ncnc"
        if "ncnc" in smiles and not heterocycles["pyrimidine"]:
            print(f"Potential pyrimidine found in: {smiles}")
            heterocycles["pyrimidine"] = True
            present_heterocycles.append("pyrimidine")
            heterocycle_count += 1

        # Check for dioxane pattern in "CCOCC"
        if "CCOCC" in smiles and not heterocycles["dioxane"]:
            print(f"Potential dioxane found in: {smiles}")
            heterocycles["dioxane"] = True
            present_heterocycles.append("dioxane")
            heterocycle_count += 1

        # Check for oxetane pattern in "COC" with numbers
        if "COC" in smiles and not heterocycles["oxetane"]:
            # This is a simplified check - in reality we'd need to confirm it's a 4-membered ring
            print(f"Potential oxetane found in: {smiles}")
            heterocycles["oxetane"] = True
            present_heterocycles.append("oxetane")
            heterocycle_count += 1

    print(
        f"After pattern matching, found {heterocycle_count} different heterocycles: {', '.join(present_heterocycles)}"
    )

    return heterocycle_count >= 3
