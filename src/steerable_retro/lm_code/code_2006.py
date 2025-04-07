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
    This function detects if a carbonyl addition reaction (like aldol) is present in the synthesis.
    """
    carbonyl_addition_found = False

    def dfs_traverse(node, depth=0):
        nonlocal carbonyl_addition_found

        # If we already found a carbonyl addition, no need to continue traversal
        if carbonyl_addition_found:
            return

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for specific carbonyl addition reactions
                carbonyl_addition_reactions = [
                    "Aldol condensation",
                    "Michael addition",
                    "Michael addition methyl",
                    "aza-Michael addition aromatic",
                    "aza-Michael addition secondary",
                    "aza-Michael addition primary",
                    "thia-Michael addition",
                    "oxa-Michael addition",
                    "Henry Reaction",
                    "Knoevenagel Condensation",
                    "Grignard_carbonyl",
                    "Grignard from aldehyde to alcohol",
                    "Grignard from ketone to alcohol",
                    "reductive amination with aldehyde",
                    "reductive amination with ketone",
                    "reductive amination with alcohol",
                    "Addition of primary amines to aldehydes/thiocarbonyls",
                    "Addition of primary amines to ketones/thiocarbonyls",
                    "Addition of secondary amines to ketones/thiocarbonyls",
                    "Addition of secondary amines to aldehydes/thiocarbonyls",
                ]

                for reaction_type in carbonyl_addition_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        carbonyl_addition_found = True
                        print(f"Found carbonyl addition reaction: {reaction_type}")
                        return

                # If specific reaction checks failed, check for general carbonyl addition pattern
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for carbonyl groups in reactants
                    carbonyl_fgs = [
                        "Aldehyde",
                        "Ketone",
                        "Ester",
                        "Carboxylic acid",
                        "Amide",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                        "Anhydride",
                    ]

                    carbonyl_found = False
                    carbonyl_reactant = None

                    for reactant in reactants:
                        for fg in carbonyl_fgs:
                            if checker.check_fg(fg, reactant):
                                carbonyl_found = True
                                carbonyl_reactant = reactant
                                print(f"Found carbonyl group ({fg}) in reactant: {reactant}")
                                break
                        if carbonyl_found:
                            break

                    # If carbonyl found in reactants, check for nucleophilic addition
                    if carbonyl_found and carbonyl_reactant:
                        # Check if this is a nucleophilic addition to carbonyl
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        product_mol = Chem.MolFromSmiles(product) if product else None
                        carbonyl_mol = (
                            Chem.MolFromSmiles(carbonyl_reactant) if carbonyl_reactant else None
                        )

                        if product_mol and carbonyl_mol:
                            # Check for typical carbonyl addition patterns in product
                            addition_patterns = [
                                Chem.MolFromSmarts("[C]-[C]-[O;H1]"),  # C-C-OH (alcohol)
                                Chem.MolFromSmarts("[C]-[C]-[N]"),  # C-C-N (amine)
                                Chem.MolFromSmarts(
                                    "[C]-[C]=[C]"
                                ),  # C-C=C (alkene from condensation)
                                Chem.MolFromSmarts("[C]-[C](=[O])-[C]"),  # C-C(=O)-C (new ketone)
                            ]

                            # Check if product has a new pattern that wasn't in the reactants
                            for pattern in addition_patterns:
                                if pattern and product_mol.HasSubstructMatch(pattern):
                                    # Check if this pattern wasn't in the reactants
                                    if not any(
                                        r and r.HasSubstructMatch(pattern) for r in reactant_mols
                                    ):
                                        carbonyl_addition_found = True
                                        print(
                                            f"Found general carbonyl addition pattern: new structural motif in product"
                                        )
                                        return

                            # Check if the carbonyl carbon has a new bond in the product
                            # This would require atom mapping analysis which is complex
                            # For simplicity, we'll check if the carbonyl group is consumed
                            for fg in carbonyl_fgs:
                                if checker.check_fg(fg, carbonyl_reactant) and not checker.check_fg(
                                    fg, product
                                ):
                                    print(f"Carbonyl group {fg} was consumed in the reaction")
                                    carbonyl_addition_found = True
                                    return

                except Exception as e:
                    print(f"Error in general carbonyl addition check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: carbonyl_addition_found = {carbonyl_addition_found}")

    return carbonyl_addition_found
