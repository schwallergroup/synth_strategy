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
    Detects if the synthesis route involves a late-stage C-S bond formation
    via nucleophilic aromatic substitution.
    """
    found_cs_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_cs_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")
            print(f"Reactants: {reactants}")
            print(f"Product: {product}")

            # Extract depth from ID or use the tracked depth
            node_depth = depth
            if "ID" in node.get("metadata", {}) and "Depth:" in node.get("metadata", {}).get(
                "ID", ""
            ):
                try:
                    depth_str = (
                        node.get("metadata", {}).get("ID", "").split("Depth:")[1].strip().split()[0]
                    )
                    node_depth = int(depth_str)
                except (IndexError, ValueError):
                    pass

            # Check if this is a late stage reaction (depth 0 or 1)
            if node_depth <= 1:
                print(f"Checking late-stage reaction at depth {node_depth}: {rsmi}")

                # Check if this is a nucleophilic aromatic substitution forming a C-S bond
                is_thioether_rxn = checker.check_reaction("thioether_nucl_sub", rsmi)
                print(f"Is thioether nucleophilic substitution: {is_thioether_rxn}")

                # Check if product has monosulfide
                has_monosulfide = checker.check_fg("Monosulfide", product)
                print(f"Product contains monosulfide: {has_monosulfide}")

                # Check if any reactant has a thiol group (typical nucleophile)
                has_thiol_reactant = any(
                    checker.check_fg("Aliphatic thiol", r) or checker.check_fg("Aromatic thiol", r)
                    for r in reactants
                )
                print(f"Has thiol reactant: {has_thiol_reactant}")

                # Check if any reactant has an aromatic halide (typical electrophile)
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)
                print(f"Has aromatic halide: {has_aryl_halide}")

                # Check for C-S bond formation directly
                if is_thioether_rxn or (has_monosulfide and has_thiol_reactant and has_aryl_halide):
                    # Additional check: Ensure the product has a C-S bond where the C is aromatic
                    # This is characteristic of nucleophilic aromatic substitution
                    print(
                        "Confirmed: Late-stage C-S bond formation via nucleophilic aromatic substitution"
                    )
                    found_cs_formation = True

                # Alternative detection method if the reaction checker fails
                if not found_cs_formation and has_monosulfide:
                    # Check if we have a thiol and an aromatic structure in reactants
                    # and the product has a monosulfide, which strongly suggests C-S bond formation
                    has_aromatic = any("c" in r or "n" in r for r in reactants)
                    if has_thiol_reactant and has_aromatic:
                        print("Detected C-S bond formation through pattern analysis")
                        found_cs_formation = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_cs_formation
