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
    Detects a synthesis involving heterocyclic compounds (thiophene and nitrogen heterocycles like imidazole, pyrimidine, etc.).
    """
    has_thiophene = False
    has_nitrogen_heterocycle = False

    # List of nitrogen-containing heterocycles to check
    nitrogen_heterocycles = [
        "imidazole",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "pyridine",
        "pyridazine",
        "purine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_thiophene, has_nitrogen_heterocycle

        # Check molecules for heterocycles
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check for thiophene
            if not has_thiophene and checker.check_ring("thiophene", mol_smiles):
                has_thiophene = True
                print(f"Detected thiophene in molecule at depth {depth}: {mol_smiles}")

            # Check for nitrogen heterocycles
            if not has_nitrogen_heterocycle:
                for heterocycle in nitrogen_heterocycles:
                    if checker.check_ring(heterocycle, mol_smiles):
                        has_nitrogen_heterocycle = True
                        print(
                            f"Detected {heterocycle} in molecule at depth {depth}: {mol_smiles}"
                        )
                        break

        # Check reactions that might form heterocycles
        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if the reaction product contains thiophene
            if not has_thiophene and checker.check_ring("thiophene", product):
                has_thiophene = True
                print(f"Detected thiophene formation in reaction at depth {depth}")

            # Check if the reaction product contains nitrogen heterocycles
            if not has_nitrogen_heterocycle:
                for heterocycle in nitrogen_heterocycles:
                    if checker.check_ring(heterocycle, product):
                        has_nitrogen_heterocycle = True
                        print(
                            f"Detected {heterocycle} formation in reaction at depth {depth}"
                        )
                        break

                # Also check for specific reactions that form nitrogen heterocycles
                nitrogen_heterocycle_reactions = [
                    "imidazole",
                    "triaryl-imidazole",
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "tetrazole_terminal",
                    "tetrazole_connect_regioisomere_1",
                    "tetrazole_connect_regioisomere_2",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                ]

                for rxn_type in nitrogen_heterocycle_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_nitrogen_heterocycle = True
                        print(f"Detected {rxn_type} reaction at depth {depth}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both heterocycles are involved
    strategy_present = has_thiophene and has_nitrogen_heterocycle
    print(f"Heterocycle-containing synthesis strategy detected: {strategy_present}")
    print(
        f"Found thiophene: {has_thiophene}, Found nitrogen heterocycle: {has_nitrogen_heterocycle}"
    )

    return strategy_present
