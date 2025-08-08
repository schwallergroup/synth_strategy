#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects if the synthesis involves early formation of a heterocycle
    (specifically at depth >= 2, which corresponds to early in the synthesis).
    """
    early_heterocycle_formation = False

    # List of heterocycles to check
    heterocycle_list = [
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
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # List of heterocycle formation reactions
    heterocycle_formation_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
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
        "3-nitrile-pyridine",
        "pyrazole",
        "oxadiazole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "imidazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal early_heterocycle_formation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Only consider early stages (depth >= 2)
                if depth >= 2:
                    # Check if this is a known heterocycle formation reaction
                    is_heterocycle_formation = False
                    for reaction_type in heterocycle_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_heterocycle_formation = True
                            print(f"Detected heterocycle formation reaction: {reaction_type}")
                            break

                    # If not a known reaction type, check for heterocycle appearance
                    if not is_heterocycle_formation:
                        # Check for heterocycle in product
                        product_heterocycles = []
                        for heterocycle in heterocycle_list:
                            has_ring = checker.check_ring(heterocycle, product_smiles)
                            if has_ring:
                                print(f"Found {heterocycle} in product")
                                product_heterocycles.append(heterocycle)

                        if product_heterocycles:
                            print(f"Product contains heterocycles: {product_heterocycles}")
                            # Check if any of these heterocycles were not in reactants
                            for heterocycle in product_heterocycles:
                                heterocycle_in_reactants = False
                                for reactant in reactants_smiles:
                                    if checker.check_ring(heterocycle, reactant):
                                        heterocycle_in_reactants = True
                                        print(
                                            f"Heterocycle {heterocycle} also found in reactant: {reactant}"
                                        )
                                        break

                                if not heterocycle_in_reactants:
                                    print(f"Early heterocycle formation detected at depth {depth}")
                                    print(f"Heterocycle: {heterocycle}")
                                    print(f"Reaction: {rsmi}")
                                    is_heterocycle_formation = True
                                    break

                    if is_heterocycle_formation:
                        early_heterocycle_formation = True
                        return  # Stop traversal once we find what we're looking for

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            if (
                not early_heterocycle_formation
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal")
    dfs_traverse(route)
    print(f"Result: {early_heterocycle_formation}")
    return early_heterocycle_formation
