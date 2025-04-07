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
    This function detects isoxazole ring formation from a β-ketoester and hydroxylamine.
    """
    isoxazole_formed = False
    hydroxylamine_used = False
    ketoester_used = False

    def dfs_traverse(node):
        nonlocal isoxazole_formed, hydroxylamine_used, ketoester_used

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for hydroxylamine or its derivatives in reactants
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() <= 10:  # Hydroxylamine and derivatives are small
                    # Check for common hydroxylamine patterns
                    hydroxylamine_patterns = ["NOH", "ONH2", "NH2O", "NHOH", "HONH2"]
                    if any(pattern in reactant for pattern in hydroxylamine_patterns):
                        print(f"Potential hydroxylamine derivative detected: {reactant}")
                        hydroxylamine_used = True
                        break

                    # Check for primary amine with oxygen nearby
                    if checker.check_fg("Primary amine", reactant) and "O" in reactant:
                        n_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
                        o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O")
                        if n_count == 1 and o_count >= 1 and mol.GetNumAtoms() <= 8:
                            print(f"Hydroxylamine-like structure detected: {reactant}")
                            hydroxylamine_used = True
                            break

            # Check for β-ketoester pattern in reactants
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    # Check for both ketone and ester functional groups
                    if checker.check_fg("Ketone", reactant) and checker.check_fg("Ester", reactant):
                        # Verify it's a β-ketoester by checking carbon chain
                        ketoester_pattern = Chem.MolFromSmarts("C(=O)CC(=O)O[C,c]")
                        alt_ketoester_pattern = Chem.MolFromSmarts(
                            "C(=O)C(C)C(=O)O[C,c]"
                        )  # For substituted versions
                        alt_ketoester_pattern2 = Chem.MolFromSmarts(
                            "C(=O)C[C,c]C(=O)O[C,c]"
                        )  # Another substituted version
                        if (
                            mol.HasSubstructMatch(ketoester_pattern)
                            or mol.HasSubstructMatch(alt_ketoester_pattern)
                            or mol.HasSubstructMatch(alt_ketoester_pattern2)
                        ):
                            ketoester_used = True
                            print(f"β-ketoester detected: {reactant}")

            # Check for isoxazole formation
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check if isoxazole is in product
                if checker.check_ring("isoxazole", product_smiles):
                    # Check if isoxazole wasn't in reactants
                    isoxazole_in_reactants = False
                    for r_smiles in reactants_smiles:
                        if checker.check_ring("isoxazole", r_smiles):
                            isoxazole_in_reactants = True
                            break

                    if not isoxazole_in_reactants:
                        print(f"Isoxazole ring formed in product: {product_smiles}")
                        isoxazole_formed = True

                        # Check for relevant reaction types that could form isoxazoles
                        cycloaddition_reactions = [
                            "Huisgen 1,3 dipolar cycloaddition",
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                        ]

                        for rxn_type in cycloaddition_reactions:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Detected cycloaddition reaction: {rxn_type}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if all conditions are met
    result = isoxazole_formed and (hydroxylamine_used or ketoester_used)
    print(f"Isoxazole formation from ketoester detected: {result}")
    print(
        f"Isoxazole formed: {isoxazole_formed}, Hydroxylamine used: {hydroxylamine_used}, Ketoester used: {ketoester_used}"
    )
    return result
