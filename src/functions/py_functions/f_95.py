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
    This function detects a convergent synthesis strategy involving heterocycle formation.
    """
    # List of heterocycles to check
    heterocycles = [
        "tetrazole",
        "triazole",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "benzotriazole",
        "furan",
        "thiophene",
        "pyrrole",
        "pyridine",
    ]

    # Heterocycle formation reaction types
    heterocycle_rxn_types = [
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "{tetrazole_terminal}",
        "{tetrazole_connect_regioisomere_1}",
        "{tetrazole_connect_regioisomere_2}",
        "{Huisgen_Cu-catalyzed_1,4-subst}",
        "{Huisgen_Ru-catalyzed_1,5_subst}",
        "{1,2,4-triazole_acetohydrazide}",
        "{1,2,4-triazole_carboxylic-acid/ester}",
        "{pyrazole}",
        "{benzimidazole_derivatives_carboxylic-acid/ester}",
        "{benzimidazole_derivatives_aldehyde}",
        "{benzothiazole}",
        "{benzoxazole_arom-aldehyde}",
        "{benzoxazole_carboxylic-acid}",
        "{thiazole}",
        "benzimidazole formation from aldehyde",
        "benzimidazole formation from acyl halide",
        "benzimidazole formation from ester/carboxylic acid",
        "benzoxazole formation from aldehyde",
        "benzoxazole formation from acyl halide",
        "benzoxazole formation from ester/carboxylic acid",
        "benzoxazole formation (intramolecular)",
        "benzothiazole formation from aldehyde",
        "benzothiazole formation from acyl halide",
        "benzothiazole formation from ester/carboxylic acid",
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "{Paal-Knorr pyrrole}",
        "{Fischer indole}",
        "{indole}",
        "{oxadiazole}",
        "A3 coupling to imidazoles",
        "Alkyne-imine cycloaddition",
    ]

    convergent_steps = []
    heterocycle_formations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for convergent step (multiple reactants)
                if len(reactants) >= 2:
                    print(
                        f"Found convergent step with {len(reactants)} fragments at depth {depth}"
                    )
                    # Track all convergent steps
                    convergent_steps.append(depth)

                    # Check if this reaction forms a heterocycle
                    heterocycle_formed = False

                    # Check if the reaction is a known heterocycle formation reaction
                    for rxn_type in heterocycle_rxn_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(
                                f"Detected heterocycle formation reaction: {rxn_type}"
                            )
                            heterocycle_formed = True
                            heterocycle_formations.append(depth)
                            break

                    # If not a known reaction type, check if a heterocycle appears in the product but not in ALL reactants
                    if not heterocycle_formed:
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, product):
                                # Check if the heterocycle is present in all reactants
                                reactants_with_heterocycle = 0
                                for reactant in reactants:
                                    if checker.check_ring(heterocycle, reactant):
                                        reactants_with_heterocycle += 1

                                # If heterocycle is in product but not in all reactants, it's likely being formed or modified
                                if reactants_with_heterocycle < len(reactants):
                                    print(
                                        f"Detected heterocycle formation or modification: {heterocycle} at depth {depth}"
                                    )
                                    heterocycle_formed = True
                                    heterocycle_formations.append(depth)
                                    break

                    # Additional check for reactions that might form heterocycles but aren't in our predefined lists
                    if not heterocycle_formed:
                        # Check if any heterocycle is in the product
                        product_heterocycles = set()
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, product):
                                product_heterocycles.add(heterocycle)

                        # Check reactants for heterocycles
                        reactant_heterocycles = set()
                        for reactant in reactants:
                            for heterocycle in heterocycles:
                                if checker.check_ring(heterocycle, reactant):
                                    reactant_heterocycles.add(heterocycle)

                        # If there's a heterocycle in the product that's not in any reactant
                        new_heterocycles = product_heterocycles - reactant_heterocycles
                        if new_heterocycles:
                            print(
                                f"Detected new heterocycle(s) in product: {new_heterocycles} at depth {depth}"
                            )
                            heterocycle_formed = True
                            heterocycle_formations.append(depth)

        # Check if this is the final target molecule and contains a heterocycle
        elif node["type"] == "mol" and depth == 0:
            for heterocycle in heterocycles:
                if checker.check_ring(heterocycle, node["smiles"]):
                    print(f"Final target molecule contains heterocycle: {heterocycle}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Strategy is present if there's at least one convergent step that forms a heterocycle in late-stage synthesis
    late_stage_convergent_heterocycle = any(
        depth in convergent_steps and depth <= 4 for depth in heterocycle_formations
    )

    print(f"Convergent steps: {convergent_steps}")
    print(f"Heterocycle formations: {heterocycle_formations}")
    print(f"Result: {late_stage_convergent_heterocycle}")
    return late_stage_convergent_heterocycle
