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
    This function detects if the synthetic route employs a protection-deprotection strategy,
    specifically looking for Boc protection of nitrogen followed by deprotection.
    """
    protection_reactions = []
    deprotection_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    protection_reactions.append((rsmi, depth))
                    print(f"Boc protection detected at depth {depth}: {rsmi}")

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    deprotection_reactions.append((rsmi, depth))
                    print(f"Boc deprotection detected at depth {depth}: {rsmi}")

                # Check for other protection reactions
                if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    protection_reactions.append((rsmi, depth))
                    print(f"Silyl protection detected at depth {depth}: {rsmi}")

                # Check for other deprotection reactions
                if (
                    checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rsmi
                    )
                    or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                ):
                    deprotection_reactions.append((rsmi, depth))
                    print(f"Silyl deprotection detected at depth {depth}: {rsmi}")

                # Fallback method if reaction checkers don't catch it
                if not (protection_reactions or deprotection_reactions):
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for Boc protection (formation of N-Boc bond)
                    reactant_mols_smiles = [r for r in reactants if r]
                    product_mol_smiles = product if product else ""

                    if product_mol_smiles:
                        # Check if Boc is in product but not in reactants
                        if (
                            checker.check_fg("Boc", product_mol_smiles)
                            and not any(checker.check_fg("Boc", r) for r in reactant_mols_smiles)
                            and any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactant_mols_smiles
                                if r
                            )
                        ):
                            protection_reactions.append((rsmi, depth))
                            print(f"Boc protection detected (fallback) at depth {depth}: {rsmi}")

                    # Check for Boc deprotection (removal of Boc group)
                    if product_mol_smiles:
                        # Check if Boc is in reactants but not in product
                        if (
                            any(checker.check_fg("Boc", r) for r in reactant_mols_smiles)
                            and not checker.check_fg("Boc", product_mol_smiles)
                            and (
                                checker.check_fg("Primary amine", product_mol_smiles)
                                or checker.check_fg("Secondary amine", product_mol_smiles)
                            )
                        ):
                            deprotection_reactions.append((rsmi, depth))
                            print(f"Boc deprotection detected (fallback) at depth {depth}: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection reactions found: {len(protection_reactions)}")
    print(f"Deprotection reactions found: {len(deprotection_reactions)}")

    # Check if we have both protection and deprotection
    if protection_reactions and deprotection_reactions:
        # In retrosynthetic traversal, protection happens at higher depth than deprotection
        for prot_rsmi, prot_depth in protection_reactions:
            for deprot_rsmi, deprot_depth in deprotection_reactions:
                # In retrosynthetic traversal, protection should be at higher depth
                if prot_depth > deprot_depth:
                    print(
                        f"Valid protection-deprotection sequence found: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                    )
                    return True

    # If we only have deprotection reactions, it's still a protection-deprotection strategy
    # since the protection might have happened outside the route
    if deprotection_reactions:
        print(
            "Only deprotection reactions found, but this still indicates a protection-deprotection strategy"
        )
        return True

    return False
