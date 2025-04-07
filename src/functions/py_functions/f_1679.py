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
    This function detects if the synthetic route involves late-stage halogenation
    (introduction of halogen atoms in the final 1-2 steps of synthesis).
    """
    halogen_introduced = False
    max_depth = 0
    final_product_smiles = route["smiles"]
    final_product = Chem.MolFromSmiles(final_product_smiles)

    # List of halogenation reaction types to check
    halogenation_reactions = [
        "Aromatic fluorination",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Chlorination",
        "Fluorination",
        "Iodination",
        "Bromination",
        "Halodeboronation of boronic acids",
        "Halodeboronation of boronic esters",
        "Wohl-Ziegler bromination benzyl primary",
        "Wohl-Ziegler bromination benzyl secondary",
        "Wohl-Ziegler bromination benzyl tertiary",
        "Wohl-Ziegler bromination allyl primary",
        "Wohl-Ziegler bromination allyl secondary",
        "Wohl-Ziegler bromination allyl tertiary",
        "Wohl-Ziegler bromination carbonyl primary",
        "Wohl-Ziegler bromination carbonyl secondary",
        "Wohl-Ziegler bromination carbonyl tertiary",
    ]

    # Additional halogen-containing functional groups
    halogen_fgs = [
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide",
        "Haloalkyne",
        "Trifluoro group",
        "Trichloro group",
    ]

    # Check if final product contains halogens
    halogen_pattern = "[F,Cl,Br,I]"
    halogen_query = Chem.MolFromSmarts(halogen_pattern)
    final_product_has_halogens = False

    if final_product:
        final_product_has_halogens = (
            len(final_product.GetSubstructMatches(halogen_query)) > 0
        )
        print(f"Final product contains halogens: {final_product_has_halogens}")

    def dfs_traverse(node, depth=0):
        nonlocal halogen_introduced, max_depth

        max_depth = max(max_depth, depth)

        # Check for late-stage halogenation (final 1-2 steps)
        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a known halogenation reaction
                for reaction_type in halogenation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Halogenation reaction detected at depth {depth}: {reaction_type}"
                        )
                        halogen_introduced = True
                        break

                # If no specific reaction type matched, check for halogen increase
                if not halogen_introduced:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Handle multiple reactants
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r
                    ]
                    product_mol = Chem.MolFromSmiles(product_part)

                    if all(reactant_mols) and product_mol:
                        # Count halogens in reactants and product
                        reactant_halogens = sum(
                            len(mol.GetSubstructMatches(halogen_query))
                            for mol in reactant_mols
                        )
                        product_halogens = len(
                            product_mol.GetSubstructMatches(halogen_query)
                        )

                        print(
                            f"Halogen count: reactants={reactant_halogens}, product={product_halogens}"
                        )

                        if product_halogens > reactant_halogens:
                            print(
                                f"Halogen increase detected at depth {depth}: {reactant_halogens} → {product_halogens}"
                            )
                            halogen_introduced = True

                        # Check for halogen-containing functional groups in product
                        if not halogen_introduced:
                            for fg in halogen_fgs:
                                if checker.check_fg(fg, product_part):
                                    # Check if this FG was not in reactants
                                    if not any(
                                        checker.check_fg(fg, r)
                                        for r in reactants_part.split(".")
                                    ):
                                        print(
                                            f"Halogen-containing functional group {fg} introduced at depth {depth}"
                                        )
                                        halogen_introduced = True
                                        break

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start DFS traversal
    dfs_traverse(route)

    # Only consider it late-stage halogenation if the final product has halogens
    result = halogen_introduced and final_product_has_halogens
    print(f"Max depth of synthesis: {max_depth}")
    print(f"Late-stage halogenation detected: {result}")
    return result
