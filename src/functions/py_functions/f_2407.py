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
    This function detects if a synthetic route involves formation of heterocycles.
    """
    # List of heterocycles to check
    heterocycles = [
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
        "isoxazole",
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

    # Track heterocycle formations
    heterocycle_formations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}")
                print(f"Product: {product_smiles}")
                print(f"Reactants: {reactants_smiles}")

                # Check which heterocycles are in the product
                product_heterocycles = []
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.append(heterocycle)
                        print(f"Product contains {heterocycle}")

                # If no heterocycles in product, skip further checks
                if not product_heterocycles:
                    print("No heterocycles found in product")
                else:
                    # Check which heterocycles are in the reactants
                    reactant_heterocycles = set()
                    for reactant_smiles in reactants_smiles:
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant_smiles):
                                reactant_heterocycles.add(heterocycle)
                                print(f"Reactant contains {heterocycle}")

                    # Find heterocycles that are in product but not in reactants
                    for heterocycle in product_heterocycles:
                        if heterocycle not in reactant_heterocycles:
                            heterocycle_formations.append((heterocycle, depth))
                            print(
                                f"Heterocycle formation detected: {heterocycle} at depth {depth}"
                            )

                            # Check if this is a known heterocycle formation reaction
                            if checker.check_reaction(
                                "Formation of NOS Heterocycles", rsmi
                            ):
                                print("Confirmed as NOS Heterocycle formation reaction")
                            elif checker.check_reaction(
                                "Paal-Knorr pyrrole synthesis", rsmi
                            ):
                                print("Confirmed as Paal-Knorr pyrrole synthesis")
                            elif checker.check_reaction(
                                "{benzimidazole_derivatives_carboxylic-acid/ester}",
                                rsmi,
                            ):
                                print(
                                    "Confirmed as benzimidazole formation from carboxylic acid/ester"
                                )
                            elif checker.check_reaction(
                                "{benzimidazole_derivatives_aldehyde}", rsmi
                            ):
                                print(
                                    "Confirmed as benzimidazole formation from aldehyde"
                                )
                            elif checker.check_reaction("{benzothiazole}", rsmi):
                                print("Confirmed as benzothiazole formation")
                            elif checker.check_reaction(
                                "{benzoxazole_arom-aldehyde}", rsmi
                            ):
                                print(
                                    "Confirmed as benzoxazole formation from aromatic aldehyde"
                                )
                            elif checker.check_reaction(
                                "{benzoxazole_carboxylic-acid}", rsmi
                            ):
                                print(
                                    "Confirmed as benzoxazole formation from carboxylic acid"
                                )
                            elif checker.check_reaction("{thiazole}", rsmi):
                                print("Confirmed as thiazole formation")
                            elif checker.check_reaction("{tetrazole_terminal}", rsmi):
                                print("Confirmed as tetrazole formation")
                            elif checker.check_reaction(
                                "{1,2,4-triazole_acetohydrazide}", rsmi
                            ):
                                print(
                                    "Confirmed as 1,2,4-triazole formation from acetohydrazide"
                                )
                            elif checker.check_reaction(
                                "{1,2,4-triazole_carboxylic-acid/ester}", rsmi
                            ):
                                print(
                                    "Confirmed as 1,2,4-triazole formation from carboxylic acid/ester"
                                )
                            elif checker.check_reaction("{pyrazole}", rsmi):
                                print("Confirmed as pyrazole formation")
                            elif checker.check_reaction("{Paal-Knorr pyrrole}", rsmi):
                                print("Confirmed as Paal-Knorr pyrrole formation")
                            elif checker.check_reaction("{Fischer indole}", rsmi):
                                print("Confirmed as Fischer indole synthesis")
                            elif checker.check_reaction("{benzofuran}", rsmi):
                                print("Confirmed as benzofuran formation")
                            elif checker.check_reaction("{benzothiophene}", rsmi):
                                print("Confirmed as benzothiophene formation")
                            elif checker.check_reaction("{indole}", rsmi):
                                print("Confirmed as indole formation")
                            elif checker.check_reaction("{oxadiazole}", rsmi):
                                print("Confirmed as oxadiazole formation")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

            # Continue traversal
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)
        else:  # node["type"] == "mol"
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found heterocycle formations
    if heterocycle_formations:
        print(f"Found heterocycle formations: {heterocycle_formations}")
        return True

    return False
