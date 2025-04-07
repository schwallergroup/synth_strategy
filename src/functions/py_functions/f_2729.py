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
    Detects if the route involves both addition and removal of nitro groups
    at different stages of the synthesis.
    """
    nitro_addition = False
    nitro_removal = False

    # Track molecules with nitro groups at different depths
    molecules_with_nitro = {}

    def dfs_traverse(node, depth=0):
        nonlocal nitro_addition, nitro_removal

        # Check molecule nodes for nitro groups
        if node["type"] == "mol":
            if node["smiles"].strip():
                try:
                    if checker.check_fg("Nitro group", node["smiles"]):
                        molecules_with_nitro[depth] = node["smiles"]
                        print(
                            f"Found molecule with nitro group at depth {depth}: {node['smiles']}"
                        )
                except Exception as e:
                    print(f"Error checking molecule for nitro group: {e}")

        # Check reaction nodes for nitro addition/removal
        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for specific nitration reactions
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                ):
                    nitro_addition = True
                    print(f"Detected nitro group addition at depth {depth}: {rsmi}")

                # Check for nitro group reduction reactions
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_removal = True
                    print(f"Detected nitro group removal at depth {depth}: {rsmi}")

                # Fallback to counting nitro groups if specific reactions not detected
                if not (nitro_addition and nitro_removal):
                    # Count nitro groups in reactants
                    reactant_has_nitro = False
                    for reactant in reactants:
                        if reactant.strip():
                            try:
                                if checker.check_fg("Nitro group", reactant):
                                    reactant_has_nitro = True
                                    break
                            except Exception as e:
                                print(f"Error checking reactant: {e}")
                                continue

                    # Check if product has nitro group
                    product_has_nitro = False
                    try:
                        if product.strip():
                            product_has_nitro = checker.check_fg("Nitro group", product)
                    except Exception as e:
                        print(f"Error checking product: {e}")

                    # In retrosynthesis, product → reactants, but we interpret in forward direction:
                    # If product has nitro but reactants don't, it's nitro addition in forward direction
                    if product_has_nitro and not reactant_has_nitro:
                        nitro_addition = True
                        print(
                            f"Detected nitro group addition at depth {depth} (by counting): {rsmi}"
                        )

                    # If reactants have nitro but product doesn't, it's nitro removal in forward direction
                    elif reactant_has_nitro and not product_has_nitro:
                        nitro_removal = True
                        print(
                            f"Detected nitro group removal at depth {depth} (by counting): {rsmi}"
                        )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have molecules with nitro groups at different depths
    # This is an additional check in case the reaction detection missed something
    if len(molecules_with_nitro) >= 2 and not (nitro_addition and nitro_removal):
        depths = sorted(molecules_with_nitro.keys())
        if max(depths) - min(depths) >= 2:  # Ensure significant depth difference
            print(f"Found molecules with nitro groups at different depths: {depths}")
            # If we have a molecule with nitro at early stage but not at late stage,
            # it suggests nitro removal
            if min(depths) <= 1 and not nitro_removal:
                nitro_removal = True
                print(
                    "Inferred nitro removal from molecule presence at different depths"
                )
            # If we have a molecule with nitro at late stage but not at early stage,
            # it suggests nitro addition
            if max(depths) >= 2 and not nitro_addition:
                nitro_addition = True
                print(
                    "Inferred nitro addition from molecule presence at different depths"
                )

    print(f"Nitro addition: {nitro_addition}, Nitro removal: {nitro_removal}")
    return nitro_addition and nitro_removal
