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
    This function detects if the synthesis involves construction of a complex
    heterocyclic system through ring formation.
    """
    heterocycle_formation_found = False

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
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if heterocycles are present in reactants
            reactant_heterocycles = set()
            for reactant in reactants:
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, reactant):
                        reactant_heterocycles.add(heterocycle)

            # Check if heterocycles are present in product
            product_heterocycles = set()
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, product):
                    product_heterocycles.add(heterocycle)

            # Check if new heterocycles are formed
            new_heterocycles = product_heterocycles - reactant_heterocycles
            if new_heterocycles:
                print(f"Heterocycle formation detected at depth {depth}: {new_heterocycles}")
                heterocycle_formation_found = True

            # Also check for specific heterocycle-forming reactions
            heterocycle_reactions = [
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
                "Fischer indole",
                "Friedlaender chinoline",
                "benzofuran",
                "benzothiophene",
                "indole",
                "oxadiazole",
                "imidazole",
            ]

            for reaction_type in heterocycle_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(
                        f"Heterocycle-forming reaction detected at depth {depth}: {reaction_type}"
                    )
                    heterocycle_formation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return heterocycle_formation_found
