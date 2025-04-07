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
    This function detects a strategy involving fluorinated building blocks
    in heterocycle synthesis.
    """
    # Track key features
    has_heterocycle_formation = False
    fluorinated_reactants = 0
    has_fluorinated_heterocycle_product = False

    # List of heterocyclic rings to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    # List of heterocycle formation reactions to check
    heterocycle_reactions = [
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
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Formation of NOS Heterocycles",
    ]

    def has_fluorine(smiles):
        """Check if a molecule contains fluorine atoms"""
        return (
            checker.check_fg("Trifluoro group", smiles)
            or (checker.check_fg("Primary halide", smiles) and "F" in smiles)
            or (checker.check_fg("Secondary halide", smiles) and "F" in smiles)
            or (checker.check_fg("Tertiary halide", smiles) and "F" in smiles)
            or (checker.check_fg("Aromatic halide", smiles) and "F" in smiles)
            or "F" in smiles
        )

    def has_heterocycle_struct(smiles):
        """Check if a molecule contains a heterocyclic structure"""
        return any(checker.check_ring(ring, smiles) for ring in heterocycles)

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation, fluorinated_reactants, has_fluorinated_heterocycle_product

        if node["type"] == "mol":
            if node["smiles"]:
                # Check if this molecule has both a heterocycle and fluorine
                if has_heterocycle_struct(node["smiles"]) and has_fluorine(node["smiles"]):
                    print(f"Found molecule with heterocycle and fluorine: {node['smiles']}")
                    # If this is the final product (depth 0), mark it
                    if depth == 0:
                        has_fluorinated_heterocycle_product = True
                        print("Final product contains fluorinated heterocycle")

                # Count fluorinated reactants/intermediates
                if has_fluorine(node["smiles"]) and depth > 0:
                    fluorinated_reactants += 1
                    print(f"Found fluorinated reactant/intermediate: {node['smiles']}")

        elif node["type"] == "reaction":
            # Check if this is a heterocycle formation reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]

                # Check for specific heterocycle formation reactions
                for reaction in heterocycle_reactions:
                    if checker.check_reaction(reaction, rxn_smiles):
                        has_heterocycle_formation = True
                        print(f"Found heterocycle formation reaction: {reaction}")
                        break

                # If no specific reaction was found, check if the reaction creates a heterocycle
                if not has_heterocycle_formation:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    # Check if product has heterocycle but reactants don't
                    if has_heterocycle_struct(product) and not all(
                        has_heterocycle_struct(r) for r in reactants
                    ):
                        has_heterocycle_formation = True
                        print(f"Found general heterocycle formation reaction: {rxn_smiles}")

                    # Check if reaction modifies a heterocycle and involves fluorine
                    if has_heterocycle_struct(product) and any(has_fluorine(r) for r in reactants):
                        print("Reaction modifies heterocycle with fluorinated building block")

                        # If the product has fluorine that wasn't in the heterocycle reactants
                        if has_fluorine(product) and not any(
                            has_heterocycle_struct(r) and has_fluorine(r) for r in reactants
                        ):
                            has_heterocycle_formation = True
                            print("Reaction adds fluorine to heterocycle")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this strategy is present
    # Modified criteria: fluorinated reactants AND either heterocycle formation OR fluorinated heterocycle product
    strategy_present = fluorinated_reactants >= 1 and (
        has_heterocycle_formation or has_fluorinated_heterocycle_product
    )

    print(f"Fluorinated heterocycle strategy detection results:")
    print(f"  Heterocycle formation reaction present: {has_heterocycle_formation}")
    print(f"  Fluorinated reactants/intermediates: {fluorinated_reactants}")
    print(f"  Fluorinated heterocycle in final product: {has_fluorinated_heterocycle_product}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
