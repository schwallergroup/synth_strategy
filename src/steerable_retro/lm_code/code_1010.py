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
    This function detects if the synthetic route involves late-stage heterocycle formation
    (thiazole formation in the second-to-last step).
    """
    heterocycle_formed = False
    depth_of_formation = -1
    max_depth = -1

    # List of heterocycles to check for formation
    heterocycles = [
        "thiazole",
        "benzothiazole",
        "benzoxazole",
        "benzimidazole",
        "oxazole",
        "isoxazole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
    ]

    # List of known heterocycle formation reactions
    heterocycle_formation_reactions = [
        "{thiazole}",
        "{benzothiazole}",
        "{benzothiazole_derivatives_carboxylic-acid/ester}",
        "{benzothiazole_derivatives_aldehyde}",
        "{benzoxazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{pyrazole}",
        "{oxadiazole}",
        "{imidazole}",
        "thiazole",
        "benzothiazole",
        "benzoxazole",
        "benzimidazole",
        "oxazole",
        "isoxazole",
        "imidazole",
        "pyrazole",
        "triazole",
        "tetrazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, depth_of_formation, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # First check if this is a known heterocycle formation reaction
                for reaction_name in heterocycle_formation_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Detected {reaction_name} formation reaction: {rsmi}")
                        heterocycle_formed = True
                        depth_of_formation = depth
                        print(f"Heterocycle formation detected at depth {depth}")
                        break

                # If no specific reaction was detected, check for heterocycle appearance
                if not heterocycle_formed:
                    for heterocycle in heterocycles:
                        # Check if reactants contain the heterocycle
                        reactants_have_heterocycle = False
                        for r_smiles in reactants_smiles.split("."):
                            if checker.check_ring(heterocycle, r_smiles):
                                reactants_have_heterocycle = True
                                print(f"Reactant contains {heterocycle}: {r_smiles}")
                                break

                        # Check if product contains the heterocycle
                        product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)
                        if product_has_heterocycle:
                            print(f"Product contains {heterocycle}: {product_smiles}")

                        # If heterocycle is in product but not in reactants, it was formed
                        if product_has_heterocycle and not reactants_have_heterocycle:
                            print(f"{heterocycle} ring formation detected at depth {depth}")
                            heterocycle_formed = True
                            depth_of_formation = depth
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Max depth: {max_depth}, Heterocycle formation depth: {depth_of_formation}")

    # Check if heterocycle formation occurred in a late stage (within first 3 steps)
    # In a retrosynthetic tree, lower depths correspond to later stages in the forward synthesis
    return heterocycle_formed and depth_of_formation <= 3 and depth_of_formation >= 0
