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
    Detects if the synthesis involves a piperazine scaffold with protection/deprotection strategy.
    """
    # Track if we found the pattern
    found_pattern = False

    # Track the depth at which we find key transformations
    boc_protection_present = False
    cbz_protection_present = False
    deprotection_depth = None
    piperazine_present = False

    # Track molecules for analysis
    protected_molecules = []
    deprotected_molecules = []

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, boc_protection_present, cbz_protection_present, deprotection_depth, piperazine_present
        nonlocal protected_molecules, deprotected_molecules

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperazine scaffold
            if checker.check_ring("piperazine", mol_smiles):
                piperazine_present = True
                print(f"Found piperazine scaffold in molecule: {mol_smiles}")

                # Check for Boc protection
                if checker.check_fg("Boc", mol_smiles):
                    boc_protection_present = True
                    protected_molecules.append(mol_smiles)
                    print(f"Found Boc protection on piperazine in molecule: {mol_smiles}")

                # Check for Cbz (benzyloxycarbonyl) protection
                # Cbz group is represented as C(=O)OCc1ccccc1
                if "C(=O)OCc1ccccc" in mol_smiles:
                    cbz_protection_present = True
                    protected_molecules.append(mol_smiles)
                    print(f"Found Cbz protection on piperazine in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            # Check for deprotection reactions
            if "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]

                # Extract reactants and product
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check if piperazine is involved in the reaction
                reactant_has_piperazine = any(
                    checker.check_ring("piperazine", r) for r in reactants
                )
                product_has_piperazine = checker.check_ring("piperazine", product)

                if reactant_has_piperazine and product_has_piperazine:
                    # Check for Boc deprotection
                    if (
                        checker.check_reaction("Boc amine deprotection", rxn_smiles)
                        or checker.check_reaction("Boc amine deprotection of guanidine", rxn_smiles)
                        or checker.check_reaction("Boc amine deprotection to NH-NH2", rxn_smiles)
                        or checker.check_reaction("Tert-butyl deprotection of amine", rxn_smiles)
                    ):
                        deprotection_depth = depth
                        deprotected_molecules.append(product)
                        print(f"Found Boc deprotection of piperazine at depth {depth}")

                    # Check for Cbz deprotection
                    elif checker.check_reaction(
                        "Carboxyl benzyl deprotection", rxn_smiles
                    ) or checker.check_reaction("Hydroxyl benzyl deprotection", rxn_smiles):
                        deprotection_depth = depth
                        deprotected_molecules.append(product)
                        print(f"Found Cbz deprotection of piperazine at depth {depth}")

                    # Check for general deprotection pattern
                    # If a protected molecule loses a protecting group
                    elif any(
                        "C(=O)OCc1ccccc" in r and "C(=O)OCc1ccccc" not in product for r in reactants
                    ) or any(
                        checker.check_fg("Boc", r) and not checker.check_fg("Boc", product)
                        for r in reactants
                    ):
                        deprotection_depth = depth
                        deprotected_molecules.append(product)
                        print(f"Found general deprotection of piperazine at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern
    if (
        piperazine_present
        and (boc_protection_present or cbz_protection_present)
        and deprotection_depth is not None
    ):
        # Deprotection should happen at a low depth (late in synthesis)
        if deprotection_depth <= 3:  # Allow slightly deeper reactions as "late-stage"
            found_pattern = True
            print("Found piperazine scaffold with protection/deprotection strategy")
            print(
                f"Protection type: {'Boc' if boc_protection_present else ''} {'Cbz' if cbz_protection_present else ''}"
            )
            print(f"Deprotection depth: {deprotection_depth}")

    return found_pattern
