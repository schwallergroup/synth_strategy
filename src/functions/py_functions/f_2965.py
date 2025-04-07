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
    This function detects a synthetic strategy involving the elaboration of a pyridine core
    with the formation of an amide linkage to a cyclohexyl moiety.
    """
    # Initialize tracking variables
    has_pyridine_core = False
    has_amide_formation = False
    has_cyclohexyl_moiety = False
    has_connected_moieties = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyridine_core, has_amide_formation, has_cyclohexyl_moiety, has_connected_moieties

        if node["type"] == "mol":
            # Check if molecule contains pyridine core
            if checker.check_ring("pyridine", node["smiles"]):
                has_pyridine_core = True
                print(f"Depth {depth}: Detected pyridine core in {node['smiles']}")

            # Check for cyclohexyl moiety
            if checker.check_ring("cyclohexane", node["smiles"]):
                has_cyclohexyl_moiety = True
                print(f"Depth {depth}: Detected cyclohexyl moiety in {node['smiles']}")

            # Check if this molecule has both pyridine and cyclohexyl with an amide linkage
            if (
                checker.check_ring("pyridine", node["smiles"])
                and checker.check_ring("cyclohexane", node["smiles"])
                and checker.check_fg("Secondary amide", node["smiles"])
            ):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # This is a strong indication that they're connected via amide
                    has_connected_moieties = True
                    print(
                        f"Depth {depth}: Detected molecule with pyridine and cyclohexyl connected via amide: {node['smiles']}"
                    )

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Schotten-Baumann_amide",
                    "Ester with primary amine to amide",
                    "Acylation of primary amines",
                    "Carboxylic acid to amide conversion",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in amide_formation_reactions
                )

                if is_amide_formation:
                    # Check if product has secondary amide
                    has_amide_in_product = checker.check_fg(
                        "Secondary amide", product_smiles
                    )

                    # Verify that one reactant has pyridine and another has cyclohexyl,
                    # or one has both and we're forming an internal amide
                    has_pyridine_reactant = any(
                        checker.check_ring("pyridine", r) for r in reactants_smiles
                    )
                    has_cyclohexyl_reactant = any(
                        checker.check_ring("cyclohexane", r) for r in reactants_smiles
                    )

                    # Check if product has both moieties
                    has_both_in_product = checker.check_ring(
                        "pyridine", product_smiles
                    ) and checker.check_ring("cyclohexane", product_smiles)

                    if (
                        has_amide_in_product
                        and has_both_in_product
                        and (has_pyridine_reactant and has_cyclohexyl_reactant)
                    ):
                        has_amide_formation = True
                        print(
                            f"Depth {depth}: Detected amide formation between pyridine and cyclohexyl moieties"
                        )

                # Check for carboxylic acid on pyridine or amine on cyclohexyl as precursors
                for reactant in reactants_smiles:
                    if checker.check_ring("pyridine", reactant) and checker.check_fg(
                        "Carboxylic acid", reactant
                    ):
                        print(
                            f"Depth {depth}: Detected pyridine with carboxylic acid - potential precursor"
                        )
                    elif checker.check_ring(
                        "cyclohexane", reactant
                    ) and checker.check_fg("Primary amine", reactant):
                        print(
                            f"Depth {depth}: Detected cyclohexyl with amine - potential precursor"
                        )

                # Check if product has both moieties connected via amide
                if has_both_in_product and has_amide_in_product:
                    has_connected_moieties = True
                    print(
                        f"Depth {depth}: Detected product with pyridine and cyclohexyl connected via amide"
                    )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present - either we found an amide formation reaction between
    # the right components, or we found a final product with both moieties connected via amide
    strategy_present = (
        has_pyridine_core
        and has_cyclohexyl_moiety
        and (has_amide_formation or has_connected_moieties)
    )

    print(
        f"Pyridine elaboration with amide linkage strategy detection result: {strategy_present}"
    )
    print(f"Detected pyridine: {has_pyridine_core}")
    print(f"Detected cyclohexyl: {has_cyclohexyl_moiety}")
    print(f"Detected amide formation: {has_amide_formation}")
    print(f"Detected connected moieties: {has_connected_moieties}")

    return strategy_present
