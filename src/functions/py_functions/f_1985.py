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
    Detects if the synthesis route involves formation of a heterocyclic system.
    """
    heterocycle_formation = False

    # List of common heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "oxirane",
        "aziridine",
        "oxetane",
        "azetidine",
        "tetrahydrofuran",
        "pyrrolidine",
        "tetrahydropyran",
        "dioxane",
    ]

    # List of reactions commonly used for heterocycle formation
    heterocycle_forming_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Pictet-Spengler",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "oxadiazole",
        "Formation of NOS Heterocycles",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known heterocycle-forming reaction
            for reaction_type in heterocycle_forming_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(
                        f"Heterocycle formation detected: {reaction_type} reaction at depth {depth}"
                    )
                    heterocycle_formation = True
                    return

            # Check if product contains heterocycles that reactants don't
            product_heterocycles = []
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, product):
                    product_heterocycles.append(heterocycle)

            if product_heterocycles:
                # Check if any reactant already has these heterocycles
                for heterocycle in product_heterocycles:
                    has_heterocycle_in_reactants = False
                    for reactant in reactants:
                        if checker.check_ring(heterocycle, reactant):
                            has_heterocycle_in_reactants = True
                            break

                    if not has_heterocycle_in_reactants:
                        print(
                            f"Heterocycle formation detected: {heterocycle} formed at depth {depth}"
                        )
                        heterocycle_formation = True
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formation
