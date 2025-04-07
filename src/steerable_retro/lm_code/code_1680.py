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
    This function detects if the synthetic route primarily involves heterocyclic systems,
    particularly nitrogen-containing heterocycles like pyrazoles and pyrimidines.
    """
    heterocycle_reactions = 0
    total_reactions = 0

    # List of heterocycles to check
    heterocycles = [
        "pyrazole",
        "pyrimidine",
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
        "imidazole",
        "oxazole",
        "thiazole",
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

    def dfs_traverse(node):
        nonlocal heterocycle_reactions, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            heterocycle_involved = False

            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                products = products_part.split(".")

                # Check if heterocycles are formed or modified in the reaction
                reactant_heterocycles = set()
                product_heterocycles = set()

                # Check heterocycles in reactants
                for reactant in reactants:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            reactant_heterocycles.add(heterocycle)

                # Check heterocycles in products
                for product in products:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product):
                            product_heterocycles.add(heterocycle)

                # Heterocycle reaction if:
                # 1. New heterocycle formed (in products but not in reactants)
                # 2. Heterocycle modified (different heterocycles in reactants vs products)
                # 3. Reaction specifically involves heterocycle chemistry

                if (
                    product_heterocycles - reactant_heterocycles
                    or reactant_heterocycles - product_heterocycles
                    or (reactant_heterocycles and product_heterocycles)
                ):
                    heterocycle_involved = True

                # Check for specific heterocycle-forming reactions
                heterocycle_forming_reactions = [
                    "Paal-Knorr pyrrole synthesis",
                    "Fischer indole",
                    "Friedlaender chinoline",
                    "benzofuran",
                    "benzothiophene",
                    "indole",
                    "oxadiazole",
                    "pyrazole",
                    "tetrazole_terminal",
                    "tetrazole_connect_regioisomere_1",
                    "tetrazole_connect_regioisomere_2",
                    "Huisgen_Cu-catalyzed_1,4-subst",
                    "Huisgen_Ru-catalyzed_1,5_subst",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                    "3-nitrile-pyridine",
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "benzothiazole",
                    "benzoxazole_arom-aldehyde",
                    "benzoxazole_carboxylic-acid",
                    "thiazole",
                ]

                for rxn_type in heterocycle_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        heterocycle_involved = True
                        break

                if heterocycle_involved:
                    heterocycle_reactions += 1
                    print(f"Heterocycle reaction detected: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Route is heterocycle-focused if most reactions involve heterocycles
    is_heterocycle_focused = total_reactions > 0 and heterocycle_reactions / total_reactions >= 0.7

    print(f"Heterocycle reactions: {heterocycle_reactions}/{total_reactions}")
    print(f"Heterocycle-focused synthesis: {is_heterocycle_focused}")

    return is_heterocycle_focused
