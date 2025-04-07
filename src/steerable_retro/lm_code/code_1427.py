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
    This function detects a synthetic strategy involving amine protection
    with Cbz or Boc group during the synthesis.
    """
    found_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for amine protection reactions
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                print(f"Found Boc protection reaction: {rsmi}")
                found_protection = True
                return

            # Manual check for Cbz protection (not in the reaction list)
            # Check if product has carbamate and reactant has amine
            has_amine_reactant = False
            has_carbamate_product = False

            for reactant in reactants:
                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                    "Secondary amine", reactant
                ):
                    print(f"Found amine in reactant: {reactant}")
                    has_amine_reactant = True
                    break

            if has_amine_reactant and checker.check_fg("Carbamic ester", product):
                print(f"Found carbamate in product: {product}")

                # Check if the carbamate is specifically a Cbz group (contains benzyl)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for benzyl group connected to carbamate
                    benzyl_carbamate = Chem.MolFromSmarts("C(=O)OC[c]1[cH][cH][cH][cH][cH]1")
                    if product_mol.HasSubstructMatch(benzyl_carbamate):
                        print(f"Confirmed Cbz protection: {rsmi}")
                        found_protection = True
                        return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Protected amine strategy found: {found_protection}")
    return found_protection
