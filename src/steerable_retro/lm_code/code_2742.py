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
    This function detects a synthetic strategy involving formylation in the late stage of synthesis
    (low depth in the retrosynthetic tree).
    """
    formylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal formylation_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only consider late-stage reactions (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Extract reactants and product
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a formylation reaction
                is_formylation_reaction = (
                    checker.check_reaction("Friedel-Crafts acylation", rsmi)
                    or checker.check_reaction("Acylation of olefines by aldehydes", rsmi)
                    or checker.check_reaction("Carbonylation with aryl formates", rsmi)
                )

                # Check for formyl group in product
                product_has_formyl = checker.check_fg(
                    "Aldehyde", product_smiles
                ) or checker.check_fg("Formaldehyde", product_smiles)

                if product_has_formyl:
                    print(f"Product has formyl group: {product_smiles}")

                    # Check if main reactant doesn't have formyl group
                    main_reactant_has_formyl = False
                    for reactant in reactants_smiles:
                        if reactant and Chem.MolFromSmiles(reactant):
                            # Skip small formylating agents
                            if len(reactant) > 5 and checker.check_fg("Aldehyde", reactant):
                                main_reactant_has_formyl = True
                                break

                    # Check for common formylating agents
                    formylating_agents = ["CN(C)C=O", "OC=O", "HC(=O)O", "HC(=O)Cl", "O=CH"]
                    has_formylating_agent = False

                    for r in reactants_smiles:
                        if not r:
                            continue
                        r_mol = Chem.MolFromSmiles(r)
                        if not r_mol:
                            continue

                        # Check if reactant is a small molecule with formyl group
                        if len(r) <= 5 and checker.check_fg("Aldehyde", r):
                            has_formylating_agent = True
                            print(f"Found formylating agent: {r}")
                            break

                        # Check against known formylating agents
                        for agent in formylating_agents:
                            agent_mol = Chem.MolFromSmiles(agent)
                            if agent_mol and r_mol.HasSubstructMatch(agent_mol):
                                has_formylating_agent = True
                                print(f"Found formylating agent: {r}")
                                break

                    # Formylation detected if:
                    # 1. Product has formyl group
                    # 2. Either it's a known formylation reaction OR
                    # 3. Main reactant doesn't have formyl AND a formylating agent is present
                    if is_formylation_reaction or (
                        not main_reactant_has_formyl and has_formylating_agent
                    ):
                        print(f"Late-stage formylation detected at depth {depth}")
                        formylation_detected = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return formylation_detected
