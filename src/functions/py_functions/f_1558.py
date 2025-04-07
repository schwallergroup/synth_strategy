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
    Detects if a Boc-protected amine is maintained throughout the synthesis.
    This means the Boc group is present in the final product and is not added
    or removed during the synthesis.
    """
    # Check if target molecule has Boc group
    if route["type"] == "mol":
        target_mol = route["smiles"]
        target_has_boc = checker.check_fg("Boc", target_mol)
        if not target_has_boc:
            print(f"Target molecule does not have Boc group: {target_mol}")
            return False

    # Track if Boc is maintained in all reactions
    boc_maintained = True

    def dfs_traverse(node, depth=0):
        nonlocal boc_maintained

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has Boc group
            product_has_boc = checker.check_fg("Boc", product)

            # Check if any reactant has Boc group
            reactant_has_boc = any(
                checker.check_fg("Boc", reactant) for reactant in reactants
            )

            # Check if Boc protection/deprotection reaction
            is_boc_protection = checker.check_reaction("Boc amine protection", rsmi)
            is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)

            # In retrosynthesis, a protection step means the Boc is added going backward,
            # which means it's removed going forward - this breaks the "maintained" strategy
            if is_boc_protection:
                print(f"Found Boc protection reaction at depth {depth}: {rsmi}")
                boc_maintained = False

            # In retrosynthesis, a deprotection step means the Boc is removed going backward,
            # which means it's added going forward - this breaks the "maintained" strategy
            if is_boc_deprotection:
                print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                boc_maintained = False

            # Check if Boc is maintained in this step
            if product_has_boc and not reactant_has_boc and not is_boc_protection:
                print(
                    f"Boc appears in product but not in reactants at depth {depth}: {rsmi}"
                )
                boc_maintained = False

            if not product_has_boc and reactant_has_boc and not is_boc_deprotection:
                print(
                    f"Boc appears in reactants but not in product at depth {depth}: {rsmi}"
                )
                boc_maintained = False

            if product_has_boc:
                print(f"Found Boc-protected amine in product at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Boc-protected amine maintained strategy detected: {boc_maintained}")
    return boc_maintained
