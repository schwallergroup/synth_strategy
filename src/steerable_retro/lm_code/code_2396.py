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
    This function detects synthetic routes involving the cleavage of polyethylene glycol chains
    into smaller fragments.
    """
    # Track if we found evidence of PEG chain cleavage
    has_peg_cleavage = False

    # Define minimum number of ether groups to consider as a PEG chain
    MIN_ETHER_COUNT = 3

    # PEG pattern - repeating ethylene glycol units
    PEG_PATTERN = "[OX2](-[CX4]-[CX4])-[CX4]-[CX4]-[OX2]"

    def is_peg_chain(mol_smiles):
        """Check if molecule contains a PEG-like structure with consecutive ethers"""
        try:
            # Check for minimum number of ethers
            ether_indices = checker.get_fg_atom_indices("Ether", mol_smiles)
            if not ether_indices or len(ether_indices) < MIN_ETHER_COUNT:
                return False

            # Check for consecutive ethers (typical PEG structure)
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return False

            # Look for ethylene glycol repeating units
            peg_pattern = Chem.MolFromSmarts(PEG_PATTERN)
            if mol.HasSubstructMatch(peg_pattern):
                print(f"Found PEG-like structure with {len(ether_indices)} ether groups")
                return True

            # Alternative check: look for multiple adjacent ether groups
            # This helps catch PEG variants
            ether_atoms = set()
            for indices_tuple in ether_indices:
                for atom_idx in indices_tuple[0]:
                    ether_atoms.add(atom_idx)

            # Check if ether oxygens are connected by short carbon chains
            connected_ethers = 0
            for i, indices_i in enumerate(ether_indices):
                o_i = indices_i[0][1]  # Oxygen atom index
                for j, indices_j in enumerate(ether_indices):
                    if i != j:
                        o_j = indices_j[0][1]  # Oxygen atom index
                        # Check if these oxygens are separated by 2-3 atoms (typical for PEG)
                        path = Chem.GetShortestPath(mol, o_i, o_j)
                        if 3 <= len(path) <= 4:  # Path length includes start and end atoms
                            connected_ethers += 1
                            break

            if connected_ethers >= MIN_ETHER_COUNT - 1:
                print(f"Found PEG-like structure with {connected_ethers} connected ether groups")
                return True

            return False
        except Exception as e:
            print(f"Error in is_peg_chain: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal has_peg_cleavage

        if has_peg_cleavage:
            return  # Early exit if we already found evidence

        if node["type"] == "reaction":
            # Extract reactants and products
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                parts = rsmi.split(">")
                if len(parts) < 3:
                    return

                reactants = parts[0].split(".")
                products = parts[2].split(".")

                # Forward direction: Check for PEG cleavage (1 reactant → multiple products)
                if len(reactants) == 1 and len(products) >= 2:
                    reactant_smiles = reactants[0]

                    # Check if reactant has PEG-like structure
                    if is_peg_chain(reactant_smiles):
                        # Check for known ether cleavage reactions
                        if (
                            checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                            or checker.check_reaction(
                                "Cleavage of methoxy ethers to alcohols", rsmi
                            )
                            or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                            or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        ):
                            print(f"Detected ether cleavage reaction")
                            has_peg_cleavage = True
                            return

                        # Count ether groups in products
                        ether_containing_products = 0
                        total_product_ethers = 0

                        for product in products:
                            if checker.check_fg("Ether", product):
                                ether_containing_products += 1
                                product_ether_indices = checker.get_fg_atom_indices(
                                    "Ether", product
                                )
                                total_product_ethers += len(product_ether_indices)
                                print(f"Product contains {len(product_ether_indices)} ether groups")

                        # Check if we have multiple products with ethers
                        if ether_containing_products >= 2:
                            print(f"Multiple ether-containing products detected")
                            has_peg_cleavage = True
                            return

                        # Check if total ether count decreased (indicating cleavage)
                        reactant_ether_indices = checker.get_fg_atom_indices(
                            "Ether", reactant_smiles
                        )
                        if (
                            total_product_ethers < len(reactant_ether_indices)
                            and ether_containing_products > 0
                        ):
                            print(
                                f"Ether count decreased: {len(reactant_ether_indices)} -> {total_product_ethers}"
                            )
                            has_peg_cleavage = True
                            return

                        # Check for alcohol formation (common in PEG cleavage)
                        alcohol_products = 0
                        for product in products:
                            if (
                                checker.check_fg("Primary alcohol", product)
                                or checker.check_fg("Secondary alcohol", product)
                                or checker.check_fg("Tertiary alcohol", product)
                            ):
                                alcohol_products += 1

                        if alcohol_products >= 2:
                            print(f"Multiple alcohol products from PEG chain cleavage detected")
                            has_peg_cleavage = True
                            return

                # Retrosynthetic direction: Check for PEG formation (multiple reactants → 1 product)
                # This appears as cleavage when traversing retrosynthetically
                if len(products) == 1 and len(reactants) >= 2:
                    product_smiles = products[0]

                    # Check if product has PEG-like structure
                    if is_peg_chain(product_smiles):
                        # Check for ether formation reactions in reverse
                        if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                            print(f"Detected ether formation reaction (retrosynthetic cleavage)")
                            has_peg_cleavage = True
                            return

                        # Count ether groups in reactants
                        ether_containing_reactants = 0
                        total_reactant_ethers = 0

                        for reactant in reactants:
                            if checker.check_fg("Ether", reactant):
                                ether_containing_reactants += 1
                                reactant_ether_indices = checker.get_fg_atom_indices(
                                    "Ether", reactant
                                )
                                total_reactant_ethers += len(reactant_ether_indices)

                        # Check if product has more ethers than sum of reactants (ether formation)
                        product_ether_indices = checker.get_fg_atom_indices("Ether", product_smiles)
                        if (
                            len(product_ether_indices) > total_reactant_ethers
                            and ether_containing_reactants > 0
                        ):
                            print(
                                f"Ether count increased: {total_reactant_ethers} -> {len(product_ether_indices)}"
                            )
                            has_peg_cleavage = True
                            return

                        # Check for alcohol consumption (common in PEG formation)
                        alcohol_reactants = 0
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                            ):
                                alcohol_reactants += 1

                        if alcohol_reactants >= 2 and ether_containing_reactants > 0:
                            print(f"Multiple alcohol reactants forming PEG chain detected")
                            has_peg_cleavage = True
                            return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Check if current molecule is a PEG
        elif node["type"] == "mol" and not node.get("in_stock", False):
            try:
                mol_smiles = node["smiles"]
                if is_peg_chain(mol_smiles):
                    # If we find a PEG molecule that's not a starting material,
                    # it might be involved in a cleavage strategy
                    print(f"Found intermediate PEG-like molecule: {mol_smiles}")
                    # We don't set has_peg_cleavage=True here because we need to confirm
                    # with an actual cleavage reaction
            except Exception as e:
                print(f"Error analyzing molecule: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"PEG chain cleavage strategy detected: {has_peg_cleavage}")
    return has_peg_cleavage
