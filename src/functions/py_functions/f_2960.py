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
    This function detects the use of Boc protection/deprotection sequence in the synthesis.
    """
    # Initialize tracking variables
    protected_molecules = {}  # Track molecules that have been Boc protected
    deprotected_molecules = {}  # Track molecules that have been Boc deprotected

    # Check if any starting materials already contain Boc groups
    def check_starting_materials(node):
        if node["type"] == "mol" and node.get("in_stock", False):
            if checker.check_fg("Carbamic ester", node["smiles"]):
                print(
                    f"Found starting material with Boc group: {node['smiles'][:20]}..."
                )
                protected_molecules[node["smiles"]] = float(
                    "inf"
                )  # Mark as protected at "infinite" depth

        for child in node.get("children", []):
            check_starting_materials(child)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for Boc protection reaction
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection with Boc anhydride", rsmi
                    )
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection of secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Boc amine protection of primary amine", rsmi
                    )
                ):

                    print(f"Detected Boc protection reaction at depth {depth}")
                    # Mark the product as Boc-protected
                    protected_molecules[product_smiles] = depth

                # Check for Boc deprotection reaction
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction(
                        "Boc amine deprotection of guanidine", rsmi
                    )
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):

                    print(f"Detected Boc deprotection reaction at depth {depth}")
                    # Check if any reactant was previously protected
                    for reactant in reactants_smiles:
                        if checker.check_fg("Carbamic ester", reactant):
                            # Mark the product as deprotected
                            deprotected_molecules[product_smiles] = depth

                # Alternative detection method: check for Boc group presence/absence
                for reactant in reactants_smiles:
                    if checker.check_fg("Carbamic ester", reactant):
                        if not checker.check_fg("Carbamic ester", product_smiles):
                            print(f"Detected Boc removal (FG change) at depth {depth}")
                            deprotected_molecules[product_smiles] = depth

                # Check if a non-Boc-containing molecule gets Boc protected
                if not any(
                    checker.check_fg("Carbamic ester", r) for r in reactants_smiles
                ):
                    if checker.check_fg("Carbamic ester", product_smiles):
                        # Verify it's a protection reaction by checking for primary/secondary amine in reactants
                        if any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants_smiles
                        ):
                            print(f"Detected Boc addition (FG change) at depth {depth}")
                            protected_molecules[product_smiles] = depth

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal to find starting materials with Boc groups
    check_starting_materials(route)

    # Start main traversal
    dfs_traverse(route)

    print(f"Protected molecules: {len(protected_molecules)}")
    print(f"Deprotected molecules: {len(deprotected_molecules)}")

    # Check if we have both protection and deprotection in the correct sequence
    has_sequence = False

    # Look for molecules that were protected and then deprotected
    for prot_mol, prot_depth in protected_molecules.items():
        for deprot_mol, deprot_depth in deprotected_molecules.items():
            # Deprotection should happen at a lower depth (later stage) than protection
            if deprot_depth < prot_depth:
                print(
                    f"Found protection at depth {prot_depth} and deprotection at depth {deprot_depth}"
                )
                has_sequence = True
                break

    # If we couldn't find a clear sequence but have both protection and deprotection,
    # still consider it a protection/deprotection strategy
    strategy_present = has_sequence or (
        len(protected_molecules) > 0 and len(deprotected_molecules) > 0
    )

    print(f"Protection/deprotection sequence detected: {strategy_present}")
    return strategy_present
