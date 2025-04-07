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
    Detects if the synthesis includes Boc protection of an amine group.

    This function traverses the synthesis route and checks for:
    1. Boc protection reactions (adding tert-butyloxycarbonyl to amines)
    2. Boc deprotection reactions (removing tert-butyloxycarbonyl from protected amines)
    3. Presence of Boc-protected intermediates

    Returns True if a Boc protection strategy is detected, False otherwise.
    """
    boc_protection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_detected

        # Skip further processing if we already found Boc protection
        if boc_protection_detected:
            return

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for Boc protection reactions directly
            is_boc_protection = any(
                [
                    checker.check_reaction("Boc amine protection", rsmi),
                    checker.check_reaction("Boc amine protection explicit", rsmi),
                    checker.check_reaction(
                        "Boc amine protection with Boc anhydride", rsmi
                    ),
                    checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi),
                    checker.check_reaction(
                        "Boc amine protection of secondary amine", rsmi
                    ),
                    checker.check_reaction(
                        "Boc amine protection of primary amine", rsmi
                    ),
                ]
            )

            # Check for Boc deprotection reactions directly
            is_boc_deprotection = any(
                [
                    checker.check_reaction("Boc amine deprotection", rsmi),
                    checker.check_reaction("Boc amine deprotection of guanidine", rsmi),
                    checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi),
                    checker.check_reaction("Tert-butyl deprotection of amine", rsmi),
                ]
            )

            # If we find a direct Boc protection or deprotection reaction, that's sufficient
            if is_boc_protection:
                print(f"Boc protection reaction detected: {rsmi}")
                boc_protection_detected = True
                return

            if is_boc_deprotection:
                print(f"Boc deprotection reaction detected: {rsmi}")
                boc_protection_detected = True
                return

            # Additional checks for Boc protection strategy
            # Check for amines in reactants and carbamates in products
            has_primary_amine = any(
                checker.check_fg("Primary amine", r) for r in reactants
            )
            has_secondary_amine = any(
                checker.check_fg("Secondary amine", r) for r in reactants
            )
            has_carbamate_in_product = checker.check_fg("Carbamic ester", product)

            # Look for evidence of Boc reagents in reactants
            boc_reagents = ["(BOC)2O", "BOC-Cl", "Boc2O", "Boc-Cl", "Boc anhydride"]
            has_boc_reagent = any(
                reagent.lower() in r.lower()
                for reagent in boc_reagents
                for r in reactants
            )

            # If we have amine reactants, carbamate products, and evidence of Boc chemistry
            if (has_primary_amine or has_secondary_amine) and has_carbamate_in_product:
                print(
                    f"Potential Boc protection: amine reactants and carbamate in product"
                )
                # Additional verification that it's specifically a Boc carbamate
                # This is a heuristic check for tert-butyl group connected to carbamate
                mol = Chem.MolFromSmiles(product)
                if mol:
                    for atom in mol.GetAtoms():
                        if (
                            atom.GetSymbol() == "C" and atom.GetDegree() == 4
                        ):  # Potential tert-butyl carbon
                            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                            if (
                                neighbors.count("C") >= 3
                            ):  # tert-butyl typically has 3 methyl groups
                                print(f"Found potential tert-butyl group in product")
                                boc_protection_detected = True
                                return

        elif node["type"] == "mol" and node["smiles"]:
            # Check if this molecule contains a Boc-protected amine
            mol_smiles = node["smiles"]
            if checker.check_fg("Carbamic ester", mol_smiles):
                print(f"Depth {depth}, Molecule contains carbamate group: {mol_smiles}")
                # Additional check to verify it's a Boc carbamate
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    for atom in mol.GetAtoms():
                        if (
                            atom.GetSymbol() == "C" and atom.GetDegree() == 4
                        ):  # Potential tert-butyl carbon
                            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                            if (
                                neighbors.count("C") >= 3
                            ):  # tert-butyl typically has 3 methyl groups
                                print(f"Found Boc-protected intermediate: {mol_smiles}")
                                boc_protection_detected = True
                                return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Boc protection strategy: {boc_protection_detected}")
    return boc_protection_detected
