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
    This function detects a synthetic strategy involving heterocyclic systems,
    particularly focusing on morpholine fused to pyridine.
    """
    # List of nitrogen heterocycles to check
    n_heterocycles = [
        "pyridine",
        "pyrrole",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "piperidine",
        "piperazine",
        "morpholine",
    ]

    # List of oxygen heterocycles to check
    o_heterocycles = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "morpholine",
    ]

    has_heterocycle = False
    has_morpholine_pyridine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle, has_morpholine_pyridine

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return

            try:
                # Check for morpholine and pyridine in the same molecule
                has_morpholine = checker.check_ring("morpholine", smiles)
                has_pyridine = checker.check_ring("pyridine", smiles)

                if has_morpholine and has_pyridine:
                    print(
                        f"Detected molecule with both morpholine and pyridine rings: {smiles}"
                    )
                    has_morpholine_pyridine = True
                    has_heterocycle = True

                # Check for any nitrogen heterocycle
                for ring in n_heterocycles:
                    if checker.check_ring(ring, smiles):
                        print(f"Detected nitrogen heterocycle ({ring}): {smiles}")
                        has_heterocycle = True
                        break

                # Check for any oxygen heterocycle if we haven't found a nitrogen one
                if not has_heterocycle:
                    for ring in o_heterocycles:
                        if checker.check_ring(ring, smiles):
                            print(f"Detected oxygen heterocycle ({ring}): {smiles}")
                            has_heterocycle = True
                            break
            except Exception as e:
                print(f"Error checking molecule {smiles}: {str(e)}")

        elif node["type"] == "reaction":
            try:
                # Check if this reaction forms a heterocycle
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains heterocycles not present in reactants
                    product_has_heterocycle = False
                    for ring in n_heterocycles + o_heterocycles:
                        if checker.check_ring(ring, product):
                            product_has_heterocycle = True

                            # Check if any reactant has this heterocycle
                            reactant_has_heterocycle = False
                            for reactant in reactants:
                                if checker.check_ring(ring, reactant):
                                    reactant_has_heterocycle = True
                                    break

                            if not reactant_has_heterocycle:
                                print(
                                    f"Detected heterocycle formation reaction: {rsmi}"
                                )
                                print(f"Formed heterocycle: {ring}")
                                has_heterocycle = True
                                break
            except Exception as e:
                print(f"Error checking reaction: {str(e)}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_heterocycle or has_morpholine_pyridine
