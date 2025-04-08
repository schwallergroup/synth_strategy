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
    Detects if the synthetic route uses a protected intermediate approach:
    1. Protection of a functional group
    2. Reactions on the protected intermediate
    3. Deprotection to reveal the original functional group
    """
    protection_events = []
    deprotection_events = []
    reaction_depths = set()

    # Define common protection/deprotection reaction types
    protection_reactions = [
        "Boc amine protection",
        "Alcohol protection with silyl ethers",
        "Protection of carboxylic acid",
        "Aldehyde or ketone acetalization",
        "Diol acetalization",
    ]

    deprotection_reactions = [
        "Boc amine deprotection",
        "Alcohol deprotection from silyl ethers",
        "Deprotection of carboxylic acid",
        "Cleavage of methoxy ethers to alcohols",
        "Acetal hydrolysis to aldehyde",
        "Acetal hydrolysis to diol",
        "Ketal hydrolysis to ketone",
        "Tert-butyl deprotection of amine",
    ]

    # Define common protecting groups
    protecting_groups = [
        "TMS ether protective group",
        "Silyl protective group",
        "Boc",
        "Acetal/Ketal",
    ]

    # Track molecules with protecting groups
    protected_molecules = {}

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events, reaction_depths, protected_molecules

        # Check for protecting groups in molecules
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            for pg in protecting_groups:
                if checker.check_fg(pg, mol_smiles):
                    if "in_stock" in node and node["in_stock"]:
                        # This is a starting material with a protecting group
                        protected_molecules[mol_smiles] = (depth, pg)
                        print(
                            f"Starting material with {pg} detected at depth {depth}: {mol_smiles}"
                        )

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reaction_depths.add(depth)

            # Extract reactants and product
            reactants_str = rsmi.split(">")[0]
            reactants = reactants_str.split(".")
            product = rsmi.split(">")[-1]

            # Check for protection reactions
            for prot_rxn in protection_reactions:
                if checker.check_reaction(prot_rxn, rsmi):
                    protection_events.append((depth, prot_rxn))
                    print(f"{prot_rxn} detected at depth {depth}")

            # Check for deprotection reactions
            for deprot_rxn in deprotection_reactions:
                if checker.check_reaction(deprot_rxn, rsmi):
                    deprotection_events.append((depth, deprot_rxn))
                    print(f"{deprot_rxn} detected at depth {depth}")

            # Check for protecting group appearance/disappearance
            for pg in protecting_groups:
                # Check if protecting group appears in product but not in all reactants
                if checker.check_fg(pg, product) and not all(
                    checker.check_fg(pg, r) for r in reactants
                ):
                    protection_events.append((depth, pg))
                    print(f"{pg} protection detected at depth {depth}")

                # Check if protecting group disappears from reactants to product
                if any(checker.check_fg(pg, r) for r in reactants) and not checker.check_fg(
                    pg, product
                ):
                    deprotection_events.append((depth, pg))
                    print(f"{pg} deprotection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there are both protection and deprotection events
    if protection_events and deprotection_events:
        print(f"Protection events: {protection_events}")
        print(f"Deprotection events: {deprotection_events}")

        # Check if there are matching protection/deprotection pairs
        for prot_depth, prot_type in protection_events:
            for deprot_depth, deprot_type in deprotection_events:
                # Match protection/deprotection types
                matched_pair = False

                # More flexible matching
                if (prot_type.replace("protection", "") in deprot_type) or (
                    deprot_type.replace("deprotection", "") in prot_type
                ):
                    matched_pair = True

                # Match protecting group types
                if "Boc" in prot_type and "Boc" in deprot_type:
                    matched_pair = True
                elif "silyl" in prot_type.lower() and "silyl" in deprot_type.lower():
                    matched_pair = True
                elif "Acetal" in prot_type.lower() and (
                    "Acetal" in deprot_type.lower() or "Ketal" in deprot_type.lower()
                ):
                    matched_pair = True

                # In retrosynthetic analysis, protection is at higher depth than deprotection
                if matched_pair and prot_depth > deprot_depth:
                    # Check if there are reactions between protection and deprotection
                    intermediate_reactions = False
                    for rxn_depth in reaction_depths:
                        if deprot_depth < rxn_depth < prot_depth:
                            intermediate_reactions = True
                            break

                    if intermediate_reactions:
                        print(
                            f"Protected intermediate approach detected with {prot_type} and {deprot_type}"
                        )
                        return True

    # Check for starting materials with protecting groups that get deprotected
    if protected_molecules and deprotection_events:
        for mol_smiles, (sm_depth, pg_type) in protected_molecules.items():
            for deprot_depth, deprot_type in deprotection_events:
                matched_pair = False

                if pg_type in deprot_type:
                    matched_pair = True

                if matched_pair and sm_depth > deprot_depth:
                    # Check if there are reactions between starting material and deprotection
                    intermediate_reactions = False
                    for rxn_depth in reaction_depths:
                        if deprot_depth < rxn_depth < sm_depth:
                            intermediate_reactions = True
                            break

                    if intermediate_reactions:
                        print(
                            f"Protected intermediate approach detected with starting material {pg_type} and {deprot_type}"
                        )
                        return True

    return False
