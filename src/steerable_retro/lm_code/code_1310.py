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
    Detects a synthetic strategy involving formation of a heterocycle in the middle of the synthesis.
    """
    heterocycle_formation_detected = False

    # List of heterocycles to check
    heterocycles = [
        "benzofuran",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "furan",
        "pyrrole",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "pyrazole",
        "quinoline",
        "isoquinoline",
        "purine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
    ]

    # Heterocycle-forming reaction types
    heterocycle_forming_reactions = [
        "benzofuran",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "benzothiazole",
        "benzimidazole_derivatives_aldehyde",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "indole",
        "thiazole",
        "pyrazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "Paal-Knorr pyrrole",
        "oxadiazole",
        "Fischer indole",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "Benzothiazole formation from aldehyde",
        "Benzothiazole formation from acyl halide",
        "Benzothiazole formation from ester/carboxylic acid",
        "Benzoxazole formation from aldehyde",
        "Benzoxazole formation from acyl halide",
        "Benzoxazole formation from ester/carboxylic acid",
        "Benzoxazole formation (intramolecular)",
        "Benzimidazole formation from aldehyde",
        "Benzimidazole formation from acyl halide",
        "Benzimidazole formation from ester/carboxylic acid",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
    ]

    # First determine the total depth of the route to identify the middle portion
    def get_max_depth(node, current_depth=0):
        if not node.get("children", []):
            return current_depth
        return max(get_max_depth(child, current_depth + 1) for child in node["children"])

    max_depth = get_max_depth(route)
    print(f"Maximum depth of route: {max_depth}")

    # Define middle portion as approximately the middle third of the synthesis
    mid_depth_min = max(1, max_depth // 3)
    mid_depth_max = min(max_depth - 1, 2 * max_depth // 3)
    print(f"Middle portion defined as depths {mid_depth_min} to {mid_depth_max}")

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this is a mid-route reaction
            if mid_depth_min <= depth <= mid_depth_max:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    product = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    # Check if product contains a heterocycle that reactants don't have
                    for heterocycle in heterocycles:
                        product_has_heterocycle = checker.check_ring(heterocycle, product)

                        if product_has_heterocycle:
                            print(f"Product contains {heterocycle}")

                            # Check if any reactant has the heterocycle
                            reactants_have_heterocycle = any(
                                checker.check_ring(heterocycle, reactant) for reactant in reactants
                            )

                            # If product has heterocycle but reactants don't, this is a heterocycle formation
                            if not reactants_have_heterocycle:
                                print(
                                    f"Mid-route {heterocycle} formation detected at depth {depth}"
                                )

                                # Check if this is a known heterocycle-forming reaction
                                reaction_identified = False
                                for rxn_type in heterocycle_forming_reactions:
                                    if checker.check_reaction(rxn_type, rsmi):
                                        print(f"Confirmed as {rxn_type} reaction")
                                        heterocycle_formation_detected = True
                                        reaction_identified = True
                                        break

                                # Even if we don't identify the specific reaction type, we've found a heterocycle formation
                                if not reaction_identified:
                                    print(
                                        f"Heterocycle formation detected but specific reaction type not identified"
                                    )
                                    heterocycle_formation_detected = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Heterocycle formation detected: {heterocycle_formation_detected}")
    return heterocycle_formation_detected
