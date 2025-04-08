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
    This function detects the use of protection/deprotection strategies in the synthetic route,
    specifically MOM protection of phenol and Boc protection of amine.
    """
    # Track protection and deprotection events
    mom_protected_molecules = set()
    mom_deprotected_molecules = set()
    boc_protected_molecules = set()
    boc_deprotected_molecules = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # Check for MOM protection (Phenol + MOM reagent -> MOM-protected phenol)
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi) or (
                any(checker.check_fg("Phenol", reactant) for reactant in reactants)
                and checker.check_fg("Ether", product)
                and not any(
                    checker.check_fg("Phenol", reactant)
                    for reactant in reactants
                    if checker.check_fg("Ether", reactant)
                )
            ):
                print(f"Potential MOM protection detected at depth {depth}")
                mom_protected_molecules.add(product)

            # Check for MOM deprotection (MOM-protected phenol -> Phenol)
            if checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi) or (
                any(checker.check_fg("Ether", reactant) for reactant in reactants)
                and checker.check_fg("Phenol", product)
                and not any(checker.check_fg("Phenol", reactant) for reactant in reactants)
            ):
                print(f"Potential MOM deprotection detected at depth {depth}")
                mom_deprotected_molecules.add(product)

            # Check for Boc protection (Primary amine + Boc reagent -> Boc-protected amine)
            if checker.check_reaction("Boc amine protection", rsmi) or (
                any(checker.check_fg("Primary amine", reactant) for reactant in reactants)
                and checker.check_fg("Carbamate", product)
            ):
                print(f"Boc protection detected at depth {depth}")
                boc_protected_molecules.add(product)

            # Check for Boc deprotection (Boc-protected amine -> Primary amine)
            if checker.check_reaction("Boc amine deprotection", rsmi) or (
                any(checker.check_fg("Carbamate", reactant) for reactant in reactants)
                and checker.check_fg("Primary amine", product)
            ):
                print(f"Boc deprotection detected at depth {depth}")
                boc_deprotected_molecules.add(product)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection for each strategy
    mom_protection_cycle = len(mom_protected_molecules) > 0 and len(mom_deprotected_molecules) > 0
    boc_protection_cycle = len(boc_protected_molecules) > 0 and len(boc_deprotected_molecules) > 0

    print(f"MOM protection cycle: {mom_protection_cycle}")
    print(f"Boc protection cycle: {boc_protection_cycle}")

    # Return true if either protection strategy is detected
    return mom_protection_cycle or boc_protection_cycle
