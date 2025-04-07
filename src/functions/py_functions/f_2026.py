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
    Detects heterocycle formation in early stages of synthesis.
    """
    found_heterocycle_formation = False
    max_depth = 0
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

    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "spiro-chromanone",
        "pyrazole",
        "phthalazinone",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "piperidine_indole",
        "imidazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
        "Pictet-Spengler",
    ]

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        if "metadata" in node and "depth" not in node["metadata"]:
            node["metadata"]["depth"] = current_depth

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    find_max_depth(route)
    early_stage_threshold = max(
        2, max_depth // 2
    )  # Consider first half of synthesis as early stage
    print(f"Max depth: {max_depth}, Early stage threshold: {early_stage_threshold}")

    def dfs_traverse(node):
        nonlocal found_heterocycle_formation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node["metadata"].get("depth", 0)

            # Check if this is an early stage reaction (high depth)
            if depth >= early_stage_threshold:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a known heterocycle formation reaction
                for reaction_type in heterocycle_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Found heterocycle formation reaction: {reaction_type} at depth {depth}"
                        )
                        found_heterocycle_formation = True
                        return

                # If not a known reaction, check for heterocycle formation by structure analysis
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycles in product that aren't in reactants
                product_heterocycles = set()
                reactant_heterocycles = set()

                # Identify heterocycles in product
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product):
                        product_heterocycles.add(heterocycle)

                # Identify heterocycles in reactants
                for reactant in reactants:
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)

                # Check if new heterocycles were formed
                new_heterocycles = product_heterocycles - reactant_heterocycles
                if new_heterocycles:
                    print(
                        f"Found new heterocycle formation at depth {depth}: {new_heterocycles}"
                    )
                    found_heterocycle_formation = True

        # Continue traversing
        for child in node.get("children", []):
            if (
                not found_heterocycle_formation
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return found_heterocycle_formation
