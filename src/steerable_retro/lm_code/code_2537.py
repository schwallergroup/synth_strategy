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
    This function detects if multiple different protection groups (Boc and Cbz)
    are used in the synthesis.
    """
    # Track protection groups
    boc_used = False
    cbz_used = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_used, cbz_used

        # Check molecule nodes
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]

            # Check for Boc protection group
            if checker.check_fg("Boc", mol_smiles):
                boc_used = True
                print(f"Boc group detected at depth {depth} in molecule: {mol_smiles[:20]}...")

            # Check for Cbz (carboxybenzyl) protection group
            # Cbz has a carbamic ester functional group with a benzene ring
            if checker.check_fg("Carbamic ester", mol_smiles) and checker.check_ring(
                "benzene", mol_smiles
            ):
                # Additional check to distinguish from Boc
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # If it has a carbamic ester but is not Boc, it's likely Cbz or similar
                    if not checker.check_fg("Boc", mol_smiles):
                        cbz_used = True
                        print(
                            f"Cbz group detected at depth {depth} in molecule: {mol_smiles[:20]}..."
                        )

        # Check reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for Boc protection/deprotection reactions
            if checker.check_reaction("Boc amine protection", rsmi) or checker.check_reaction(
                "Boc amine deprotection", rsmi
            ):
                boc_used = True
                print(f"Boc protection/deprotection reaction detected at depth {depth}")

            # Check for Cbz protection/deprotection reactions
            # Since there's no direct checker for Cbz reactions, look for keywords
            if "Cbz" in rsmi or "CBz" in rsmi or "benzyloxycarbonyl" in rsmi.lower():
                cbz_used = True
                print(f"Cbz protection/deprotection reaction detected at depth {depth}")

            # Check product and reactants for protection groups
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if Boc is introduced or removed in the reaction
                product_has_boc = checker.check_fg("Boc", product)
                reactants_have_boc = any(checker.check_fg("Boc", r) for r in reactants)
                if product_has_boc or reactants_have_boc:
                    boc_used = True
                    print(f"Boc group involved in reaction at depth {depth}")

                # Check if Cbz is introduced or removed in the reaction
                product_has_cbz = checker.check_fg(
                    "Carbamic ester", product
                ) and checker.check_ring("benzene", product)
                reactants_have_cbz = any(
                    checker.check_fg("Carbamic ester", r) and checker.check_ring("benzene", r)
                    for r in reactants
                )

                # Additional check to distinguish from Boc
                if (product_has_cbz and not checker.check_fg("Boc", product)) or any(
                    checker.check_fg("Carbamic ester", r)
                    and checker.check_ring("benzene", r)
                    and not checker.check_fg("Boc", r)
                    for r in reactants
                ):
                    cbz_used = True
                    print(f"Cbz group involved in reaction at depth {depth}")
            except Exception as e:
                print(f"Error analyzing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both protection groups are used
    if boc_used and cbz_used:
        print("Multiple protection groups (Boc and Cbz) strategy detected")
        return True

    if boc_used:
        print("Only Boc protection group detected")
    elif cbz_used:
        print("Only Cbz protection group detected")
    else:
        print("No protection groups detected")

    return False
