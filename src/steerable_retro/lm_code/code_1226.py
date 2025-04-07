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
    Detects if the synthesis route maintains a Boc-protected amine throughout.
    This means every reaction in the route must have Boc protection in both
    reactants and products.
    """
    reaction_count = 0
    main_pathway_molecules = []  # Track molecules in the main synthetic pathway

    def dfs_traverse(node, depth=0, is_main_pathway=True):
        nonlocal reaction_count

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Only consider molecules in the main synthetic pathway
            if is_main_pathway:
                main_pathway_molecules.append(mol_smiles)
                has_boc = checker.check_fg("Boc", mol_smiles)

                if has_boc:
                    print(f"Depth {depth}: Molecule has Boc protection: {mol_smiles}")
                else:
                    print(f"Depth {depth}: Molecule without Boc protection: {mol_smiles}")
                    # If any molecule in the main pathway doesn't have Boc, strategy fails
                    if not is_auxiliary_reagent(mol_smiles):
                        return False

            # If this is a starting material, no need to check further
            if node.get("in_stock", False):
                return True

        elif node["type"] == "reaction":
            reaction_count += 1
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth}: Checking reaction: {rsmi}")

                # Check if this is a Boc protection or deprotection reaction
                is_boc_protection = checker.check_reaction("Boc amine protection", rsmi)
                is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)

                # If this is a Boc protection/deprotection reaction in the main pathway, strategy fails
                if is_main_pathway and (is_boc_protection or is_boc_deprotection):
                    print(
                        f"Depth {depth}: Found Boc {'protection' if is_boc_protection else 'deprotection'} reaction"
                    )
                    return False

                # Check if product has Boc protection
                product_has_boc = checker.check_fg("Boc", product)

                # Check if main reactants have Boc protection (excluding auxiliary reagents)
                main_reactants = [r for r in reactants if not is_auxiliary_reagent(r)]
                reactants_have_boc = (
                    all(checker.check_fg("Boc", reactant) for reactant in main_reactants)
                    if main_reactants
                    else False
                )

                # If either the product or main reactants don't have Boc in a main pathway reaction, strategy fails
                if is_main_pathway and (not product_has_boc or not reactants_have_boc):
                    if not product_has_boc:
                        print(f"Depth {depth}: Product without Boc protection: {product}")
                    if not reactants_have_boc:
                        print(f"Depth {depth}: Not all main reactants have Boc protection")
                        for reactant in main_reactants:
                            has_boc = checker.check_fg("Boc", reactant)
                            print(f"  Reactant: {reactant}, Has Boc: {has_boc}")
                    return False

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        # For reactions, all children are part of the main pathway if the current node is
        # For molecules, only one child (the main reaction) is part of the main pathway
        if node["type"] == "reaction":
            for child in node.get("children", []):
                result = dfs_traverse(child, depth + 1, is_main_pathway)
                if not result:
                    return False
        else:  # molecule node
            for i, child in enumerate(node.get("children", [])):
                # For molecule nodes, only the first child is considered part of the main pathway
                child_is_main = is_main_pathway and i == 0
                result = dfs_traverse(child, depth + 1, child_is_main)
                if not result:
                    return False

        return True

    # Start traversal
    result = dfs_traverse(route)

    # Check if we have any reactions
    if reaction_count == 0:
        print("No reactions found in the route")
        return False

    # Check if all main pathway molecules have Boc protection
    molecules_without_boc = [
        mol
        for mol in main_pathway_molecules
        if not checker.check_fg("Boc", mol) and not is_auxiliary_reagent(mol)
    ]

    if molecules_without_boc:
        print(
            f"Found {len(molecules_without_boc)} molecules in main pathway without Boc protection"
        )
        return False

    print(f"All molecules in main pathway have Boc protection across {reaction_count} reactions")
    return True
