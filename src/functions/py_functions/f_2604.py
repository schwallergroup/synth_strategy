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
    This function detects a strategy involving late-stage functionalization
    at a benzylic position.
    """
    has_benzylic_functionalization = False
    benzylic_depth = -1

    # List of functional groups that could be involved in benzylic functionalization
    benzylic_fgs = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Ether",
        "Ester",
        "Amide",
        "Nitrile",
        "Nitro group",
        "Carboxylic acid",
        "Aldehyde",
        "Ketone",
    ]

    # List of reactions commonly used for benzylic functionalization
    benzylic_rxn_types = [
        "Wohl-Ziegler bromination benzyl primary",
        "Wohl-Ziegler bromination benzyl secondary",
        "Wohl-Ziegler bromination benzyl tertiary",
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of alcohol to carboxylic acid",
        "Reduction of ester to primary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Friedel-Crafts alkylation",
        "Friedel-Crafts acylation",
        "Aromatic chlorination",
        "Aromatic bromination",
        "Aromatic iodination",
        "Aromatic fluorination",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of ketone to carboxylic acid",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_benzylic_functionalization, benzylic_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for specific benzylic functionalization reactions first
                is_benzylic_functionalization = False
                for rxn_type in benzylic_rxn_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_benzylic_functionalization = True
                        print(
                            f"Detected benzylic functionalization reaction: {rxn_type} at depth {depth}"
                        )
                        break

                if not is_benzylic_functionalization:
                    # Convert to RDKit molecules
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(reactants):
                        # Define benzylic pattern - carbon attached to aromatic ring
                        benzylic_pattern = Chem.MolFromSmarts(
                            "[c;r6]-[C;!$(C=O);!$(C#N);!$(C=C);!$(C=N)]"
                        )

                        # Check if product has benzylic functional group
                        product_has_benzylic_fg = False

                        # First check if product has benzylic carbon
                        if product.HasSubstructMatch(benzylic_pattern):
                            benzylic_matches = product.GetSubstructMatches(
                                benzylic_pattern
                            )
                            benzylic_carbons = [
                                match[1] for match in benzylic_matches
                            ]  # The carbon in c-C

                            # Then check if any functional group is present in the product
                            for fg in benzylic_fgs:
                                if checker.check_fg(fg, product_smiles):
                                    product_has_benzylic_fg = True
                                    print(
                                        f"Product has benzylic structure and {fg} at depth {depth}"
                                    )
                                    break

                        # Check if any reactant has a benzylic carbon but different functional groups
                        reactant_has_benzylic_carbon = False
                        if product_has_benzylic_fg:
                            for r_mol, r_smiles in zip(reactants, reactants_smiles):
                                if r_mol and r_mol.HasSubstructMatch(benzylic_pattern):
                                    # Check if the reactant has different functional groups than the product
                                    reactant_has_benzylic_carbon = True

                                    # Check if the functional groups are different
                                    reactant_fgs = set()
                                    product_fgs = set()

                                    for fg in benzylic_fgs:
                                        if checker.check_fg(fg, r_smiles):
                                            reactant_fgs.add(fg)
                                        if checker.check_fg(fg, product_smiles):
                                            product_fgs.add(fg)

                                    if reactant_fgs != product_fgs:
                                        print(
                                            f"Reactant has different benzylic functional groups: {r_smiles}"
                                        )
                                        is_benzylic_functionalization = True
                                        break

                        # If we didn't find a specific pattern but both conditions are met, it's likely a functionalization
                        if (
                            product_has_benzylic_fg
                            and reactant_has_benzylic_carbon
                            and not is_benzylic_functionalization
                        ):
                            is_benzylic_functionalization = True
                            print(
                                f"Detected general benzylic functionalization at depth {depth}"
                            )

                if is_benzylic_functionalization:
                    has_benzylic_functionalization = True
                    # Update to track the latest (lowest depth) benzylic functionalization
                    if benzylic_depth == -1 or depth < benzylic_depth:
                        benzylic_depth = depth
                    print(f"Confirmed benzylic functionalization at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if benzylic functionalization is late-stage (low depth)
    # Consider depth <= 3 as late-stage based on the output
    late_stage = benzylic_depth <= 3 if benzylic_depth != -1 else False

    strategy_present = has_benzylic_functionalization and late_stage
    print(
        f"Late-stage benzylic functionalization strategy detected: {strategy_present}"
    )
    print(f"Benzylic functionalization depth: {benzylic_depth}")
    return strategy_present
