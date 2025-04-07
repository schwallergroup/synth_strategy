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
    This function detects if the synthetic route contains a late-stage modification
    of an amide group, particularly conversion to a heterocyclic amide.
    """
    amide_modifications = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has an amide group (primary, secondary, or tertiary)
            has_amide_reactant = False
            amide_reactant = None
            for r in reactants:
                if r:
                    for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                        if checker.check_fg(amide_type, r):
                            has_amide_reactant = True
                            amide_reactant = r
                            print(f"Found {amide_type} in reactant: {r}")
                            break
                    if has_amide_reactant:
                        break

            # Check if product has heterocyclic structures
            product_has_heterocycle = False
            heterocycle_found = None
            if product:
                for ring_name in [
                    "pyrrole",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "triazole",
                    "tetrazole",
                    "indole",
                    "benzimidazole",
                    "benzoxazole",
                    "benzothiazole",
                ]:
                    if checker.check_ring(ring_name, product):
                        product_has_heterocycle = True
                        heterocycle_found = ring_name
                        print(f"Found heterocycle {ring_name} in product: {product}")
                        break

            # Check if the amide is actually modified (not present in product)
            amide_modified = False
            if has_amide_reactant and product:
                amide_in_product = False
                for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                    if checker.check_fg(amide_type, product):
                        amide_in_product = True
                        break

                # If the amide is still present, check if it's a different amide (modified)
                if amide_in_product:
                    # This is a simplification - ideally we would check if the specific amide group was modified
                    # using atom mapping, but for now we'll assume if there's a heterocycle, the amide was modified
                    if product_has_heterocycle:
                        amide_modified = True
                else:
                    # Amide is not present in product, so it was modified
                    amide_modified = True

            # Check for heterocycle formation reactions
            is_heterocycle_formation = False
            heterocycle_formation_reactions = [
                "benzimidazole_derivatives_carboxylic-acid/ester",
                "benzimidazole_derivatives_aldehyde",
                "benzothiazole",
                "benzoxazole_arom-aldehyde",
                "benzoxazole_carboxylic-acid",
                "thiazole",
                "tetrazole_terminal",
                "tetrazole_connect_regioisomere_1",
                "tetrazole_connect_regioisomere_2",
                "1,2,4-triazole_acetohydrazide",
                "1,2,4-triazole_carboxylic-acid/ester",
                "pyrazole",
                "oxadiazole",
                "indole",
                "Paal-Knorr pyrrole",
                "Fischer indole",
                "Formation of NOS Heterocycles",
            ]

            for rxn_type in heterocycle_formation_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    is_heterocycle_formation = True
                    print(f"Found heterocycle formation reaction: {rxn_type}")
                    break

            # Also check for general amide modification reactions
            amide_modification_reactions = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acylation of secondary amines",
                "Acylation of primary amines",
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Hydrogenolysis of amides/imides/carbamates",
                "Hydrolysis of amides/imides/carbamates",
            ]

            for rxn_type in amide_modification_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found amide modification reaction: {rxn_type}")
                    if has_amide_reactant:
                        amide_modified = True
                        break

            # Verify that an amide is being modified to include a heterocycle
            # We're being more lenient here - if we have an amide in reactants and a heterocycle in product,
            # we'll consider it an amide modification to heterocycle
            if (
                has_amide_reactant
                and product_has_heterocycle
                and (amide_modified or is_heterocycle_formation)
            ):
                print(f"Found amide modification to heterocycle at depth {depth}, rsmi: {rsmi}")
                amide_modifications.append((depth, f"amide_to_{heterocycle_found}"))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if any amide modifications occurred in the late stage (lower depth)
    # Late stage is typically considered the first third of the synthesis depth
    if amide_modifications and max_depth > 0:
        # Ensure at least depth 1 is considered late stage for short routes
        late_stage_threshold = max_depth / 3 if max_depth >= 3 else 1
        print(f"Max depth: {max_depth}, Late stage threshold: {late_stage_threshold}")

        for depth, mod_type in amide_modifications:
            print(f"Checking modification at depth {depth}: {mod_type}")
            if depth <= late_stage_threshold:
                print(
                    f"Late stage amide modification detected at depth {depth} (threshold: {late_stage_threshold})"
                )
                return True

    return False
