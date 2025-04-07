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
    Detects if the synthesis route involves phenol alkylation to form aryl ether.
    """
    phenol_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_alkylation_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for known phenol alkylation reactions
                if (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction(
                        "Williamson Ether Synthesis (intra to epoxy)", rsmi
                    )
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Chan-Lam etherification", rsmi)
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                    )
                    or checker.check_reaction("{Williamson ether}", rsmi)
                    or checker.check_reaction(
                        "S-alkylation of thiols with alcohols", rsmi
                    )
                ):

                    print(
                        f"Found potential phenol alkylation reaction at depth {depth}"
                    )

                    # Verify that a phenol is involved in the reaction
                    phenol_in_reactants = any(
                        checker.check_fg("Phenol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )
                    ether_in_product = checker.check_fg("Ether", product)

                    if phenol_in_reactants and ether_in_product:
                        print(f"Confirmed phenol alkylation at depth {depth}")
                        phenol_alkylation_found = True

                # Check for O-alkylation reactions that might not be categorized specifically
                elif not phenol_alkylation_found:
                    # Check for phenol or aromatic alcohol in reactants
                    phenol_reactant_indices = [
                        i
                        for i, r in enumerate(reactants)
                        if checker.check_fg("Phenol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                    ]

                    # Check for alkyl halides or other alkylating agents in reactants
                    alkylating_agent_present = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Triflate", r)
                        or checker.check_fg("Mesylate", r)
                        or checker.check_fg("Tosylate", r)
                        for r in reactants
                    )

                    if phenol_reactant_indices and alkylating_agent_present:
                        print(
                            f"Found phenol and alkylating agent in reactants at depth {depth}"
                        )

                        # Check for ether in product
                        if checker.check_fg("Ether", product):
                            print(f"Found ether in product at depth {depth}")

                            # Check if product contains an aromatic ring
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                # Check if the product has an aromatic ring
                                has_aromatic = False
                                for atom in product_mol.GetAtoms():
                                    if atom.GetIsAromatic():
                                        has_aromatic = True
                                        break

                                if has_aromatic:
                                    print(
                                        f"Confirmed aryl ether formation at depth {depth}"
                                    )
                                    phenol_alkylation_found = True

                # Final check for any reaction that converts a phenol to an aryl ether
                elif not phenol_alkylation_found:
                    # Check for phenol in reactants
                    phenol_in_reactants = any(
                        checker.check_fg("Phenol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )

                    if phenol_in_reactants:
                        print(f"Found phenol in reactants at depth {depth}")

                        # Check for ether in product and no phenol in product
                        if (
                            checker.check_fg("Ether", product)
                            and not checker.check_fg("Phenol", product)
                            and not checker.check_fg("Aromatic alcohol", product)
                        ):

                            # Check if product contains an aromatic ring
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                has_aromatic = False
                                for atom in product_mol.GetAtoms():
                                    if atom.GetIsAromatic():
                                        has_aromatic = True
                                        break

                                if has_aromatic:
                                    print(
                                        f"Confirmed phenol conversion to aryl ether at depth {depth}"
                                    )
                                    phenol_alkylation_found = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return phenol_alkylation_found
