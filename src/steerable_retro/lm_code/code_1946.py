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
    This function detects heterocycle formation in the early stages of synthesis.
    Specifically looking for ring formation at high depth values (early in synthesis).
    """
    heterocycle_formed = False

    # List of common heterocycles to check
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

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains heterocycles not present in reactants
                product_heterocycles = set()
                reactants_heterocycles = set()

                # Check heterocycles in product
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_heterocycles.add(heterocycle)

                # Check heterocycles in reactants
                for reactant in reactants_smiles:
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, reactant):
                            reactants_heterocycles.add(heterocycle)

                # Check if new heterocycles were formed
                new_heterocycles = product_heterocycles - reactants_heterocycles
                if new_heterocycles:
                    print(f"Heterocycle formation detected at depth {depth}: {new_heterocycles}")

                    # Check for specific heterocycle formation reactions
                    heterocycle_formation_reactions = [
                        "Formation of NOS Heterocycles",
                        "Paal-Knorr pyrrole synthesis",
                        "benzothiazole formation from aldehyde",
                        "benzothiazole formation from acyl halide",
                        "benzothiazole formation from ester/carboxylic acid",
                        "benzoxazole formation from aldehyde",
                        "benzoxazole formation from acyl halide",
                        "benzoxazole formation from ester/carboxylic acid",
                        "benzoxazole formation (intramolecular)",
                        "benzimidazole formation from aldehyde",
                        "benzimidazole formation from acyl halide",
                        "benzimidazole formation from ester/carboxylic acid",
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        "Huisgen 1,3 dipolar cycloaddition",
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                        "Pyrazole formation",
                        "A3 coupling to imidazoles",
                        "Alkyne-imine cycloaddition",
                        "Azide-nitrile click cycloaddition to tetrazole",
                        "Azide-nitrile click cycloaddition to triazole",
                        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                        "Intramolecular amination (heterocycle formation)",
                        "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                    ]

                    for reaction_type in heterocycle_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Confirmed heterocycle formation reaction: {reaction_type}")
                            heterocycle_formed = True
                            return

                    # If no specific reaction type matched but we found new heterocycles,
                    # still mark as heterocycle formation
                    heterocycle_formed = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_formed
