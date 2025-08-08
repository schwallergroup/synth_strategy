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
    This function detects a linear synthesis route featuring SNAr C-N bond formation,
    nitro reduction, amide coupling, and late-stage heterocycle construction.
    """
    # Initialize flags for each key transformation
    has_snar = False
    has_nitro_reduction = False
    has_amide_formation = False
    has_heterocycle_formation = False
    has_methoxy = False
    has_cyclopropyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar, has_nitro_reduction, has_amide_formation, has_heterocycle_formation, has_methoxy, has_cyclopropyl

        if node["type"] == "mol" and depth == 0:
            # Check for methoxy group in final product
            if checker.check_fg("Ether", node["smiles"]):
                has_methoxy = True
                print(f"Methoxy group detected in final product: {node['smiles']}")

            # Check for cyclopropyl group in final product
            if checker.check_ring("cyclopropane", node["smiles"]):
                has_cyclopropyl = True
                print(f"Cyclopropyl group detected in final product: {node['smiles']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check for SNAr reaction - expanded to include more reaction types
            if (
                checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
            ):
                has_snar = True
                print(f"SNAr reaction detected: {rsmi}")

            # Alternative check for SNAr: look for nitro group and amine nucleophile
            if not has_snar:
                reactants = reactants_str.split(".")
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant) and any(
                        checker.check_fg(fg, reactants_str)
                        for fg in ["Primary amine", "Secondary amine", "Aniline"]
                    ):
                        has_snar = True
                        print(f"SNAr reaction detected via functional groups: {rsmi}")
                        break

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                has_nitro_reduction = True
                print(f"Nitro reduction detected: {rsmi}")

            # Alternative check for nitro reduction: nitro in reactants, amine in products
            if not has_nitro_reduction:
                if checker.check_fg("Nitro group", reactants_str) and checker.check_fg(
                    "Primary amine", product_str
                ):
                    has_nitro_reduction = True
                    print(f"Nitro reduction detected via functional groups: {rsmi}")

            # Check for amide formation - expanded list
            amide_reactions = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Ester with primary amine to amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acyl chloride with secondary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
            ]
            if any(checker.check_reaction(rxn, rsmi) for rxn in amide_reactions):
                has_amide_formation = True
                print(f"Amide formation detected: {rsmi}")

            # Alternative check for amide formation: acid/acyl halide + amine → amide
            if not has_amide_formation:
                if (
                    (
                        checker.check_fg("Carboxylic acid", reactants_str)
                        or checker.check_fg("Acyl halide", reactants_str)
                    )
                    and (
                        checker.check_fg("Primary amine", reactants_str)
                        or checker.check_fg("Secondary amine", reactants_str)
                    )
                    and checker.check_fg("Primary amide", product_str)
                    or checker.check_fg("Secondary amide", product_str)
                ):
                    has_amide_formation = True
                    print(f"Amide formation detected via functional groups: {rsmi}")

            # Check for heterocycle formation - expanded list
            heterocycle_reactions = [
                "Formation of NOS Heterocycles",
                "benzimidazole_derivatives_carboxylic-acid/ester",
                "benzimidazole_derivatives_aldehyde",
                "benzothiazole",
                "benzoxazole_arom-aldehyde",
                "benzoxazole_carboxylic-acid",
                "thiazole",
                "tetrazole_terminal",
                "1,2,4-triazole_acetohydrazide",
                "1,2,4-triazole_carboxylic-acid/ester",
                "pyrazole",
                "oxadiazole",
                "Paal-Knorr pyrrole synthesis",
                "Pictet-Spengler",
                "Fischer indole",
                "Niementowski_quinazoline",
                "Intramolecular amination (heterocycle formation)",
                "Intramolecular amination of azidobiphenyls (heterocycle formation)",
            ]
            if any(checker.check_reaction(rxn, rsmi) for rxn in heterocycle_reactions):
                has_heterocycle_formation = True
                print(f"Heterocycle formation detected: {rsmi}")

            # Alternative check for heterocycle formation: count rings in reactants vs products
            if not has_heterocycle_formation:
                # Check if a new heterocycle appears in the product
                heterocycles = [
                    "pyrrole",
                    "pyridine",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "triazole",
                    "tetrazole",
                    "indole",
                    "benzimidazole",
                    "benzoxazole",
                    "benzothiazole",
                ]

                for ring in heterocycles:
                    if not any(
                        checker.check_ring(ring, r) for r in reactants_str.split(".")
                    ) and checker.check_ring(ring, product_str):
                        has_heterocycle_formation = True
                        print(f"Heterocycle formation detected via ring check ({ring}): {rsmi}")
                        break

            # Additional check for methoxy and cyclopropyl in products
            if checker.check_fg("Ether", product_str):
                has_methoxy = True
                print(f"Methoxy group detected in product: {product_str}")

            if checker.check_ring("cyclopropane", product_str):
                has_cyclopropyl = True
                print(f"Cyclopropyl group detected in product: {product_str}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present (requires at least 3 of the key transformations)
    strategy_score = sum(
        [
            has_snar,
            has_nitro_reduction,
            has_amide_formation,
            has_heterocycle_formation,
            has_methoxy,
            has_cyclopropyl,
        ]
    )

    print(f"Strategy score: {strategy_score}/6")
    print(
        f"SNAr: {has_snar}, Nitro reduction: {has_nitro_reduction}, Amide formation: {has_amide_formation}"
    )
    print(
        f"Heterocycle formation: {has_heterocycle_formation}, Methoxy: {has_methoxy}, Cyclopropyl: {has_cyclopropyl}"
    )

    return strategy_score >= 3
