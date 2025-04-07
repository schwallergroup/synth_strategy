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
    This function detects a strategy involving olefin coupling (like Heck reaction)
    in the synthesis of polycyclic compounds.
    """
    olefin_coupling_detected = False
    polycyclic_product_detected = False

    def dfs_traverse(node):
        nonlocal olefin_coupling_detected, polycyclic_product_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for olefin coupling reactions using the checker function
                if (
                    checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                    or checker.check_reaction("Heck_non-terminal_vinyl", rsmi)
                    or checker.check_reaction("Oxidative Heck reaction", rsmi)
                    or checker.check_reaction("Oxidative Heck reaction with vinyl ester", rsmi)
                    or checker.check_reaction("Heck reaction with vinyl ester and amine", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Stille", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Negishi", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                ):
                    olefin_coupling_detected = True
                    print(f"Olefin coupling reaction detected: {rsmi}")

                # If no specific olefin coupling reaction is detected, check for general patterns
                if not olefin_coupling_detected:
                    reactants_smiles = rsmi.split(">")[0]
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if all(r is not None for r in reactants) and product is not None:
                        # Check for olefin coupling components
                        reactants_have_alkene = any(
                            checker.check_fg("Alkene", Chem.MolToSmiles(r))
                            or checker.check_fg("Vinyl", Chem.MolToSmiles(r))
                            or checker.check_fg("Allyl", Chem.MolToSmiles(r))
                            for r in reactants
                            if r is not None
                        )

                        reactants_have_coupling_partner = any(
                            checker.check_fg("Aromatic halide", Chem.MolToSmiles(r))
                            or checker.check_fg("Boronic acid", Chem.MolToSmiles(r))
                            or checker.check_fg("Boronic ester", Chem.MolToSmiles(r))
                            or checker.check_fg("Triflate", Chem.MolToSmiles(r))
                            for r in reactants
                            if r is not None
                        )

                        # Check for conjugated system in product
                        product_has_conjugated = (
                            checker.check_fg("Vinyl", product_smiles)
                            or checker.check_fg("Alkene", product_smiles)
                        ) and any(
                            checker.check_ring(ring, product_smiles)
                            for ring in [
                                "benzene",
                                "pyridine",
                                "furan",
                                "thiophene",
                                "pyrrole",
                                "naphthalene",
                                "indole",
                                "quinoline",
                                "isoquinoline",
                            ]
                        )

                        # Less restrictive condition for olefin coupling detection
                        if (
                            reactants_have_alkene and reactants_have_coupling_partner
                        ) or product_has_conjugated:
                            olefin_coupling_detected = True
                            print(f"General olefin coupling pattern detected: {rsmi}")

                # Check if product is polycyclic
                product = Chem.MolFromSmiles(product_smiles)
                if product is not None:
                    ring_info = product.GetRingInfo()
                    num_rings = ring_info.NumRings()

                    # Check for at least 2 rings (more inclusive definition of polycyclic)
                    if num_rings >= 2:
                        # Check if any rings are fused (share bonds)
                        bonds_in_rings = {}
                        for ring_idx in range(num_rings):
                            for bond_idx in ring_info.BondMembers(ring_idx):
                                if bond_idx in bonds_in_rings:
                                    # This bond is in multiple rings, indicating fused rings
                                    polycyclic_product_detected = True
                                    print(
                                        f"Polycyclic product detected with {num_rings} rings, including fused rings"
                                    )
                                    break
                                bonds_in_rings[bond_idx] = ring_idx
                            if polycyclic_product_detected:
                                break

                        # If we didn't find fused rings but have 2+ rings, still consider it polycyclic
                        if not polycyclic_product_detected:
                            polycyclic_product_detected = True
                            print(f"Polycyclic product detected with {num_rings} separate rings")

            except Exception as e:
                print(f"Error processing reaction: {str(e)}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the strategy is detected
    return olefin_coupling_detected and polycyclic_product_detected
