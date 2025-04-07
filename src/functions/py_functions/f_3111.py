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
    Detects if the synthesis involves a nitrile intermediate that is transformed
    into a heterocycle in the final steps.
    """
    has_nitrile_intermediate = False
    nitrile_to_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrile_intermediate, nitrile_to_heterocycle

        if node["type"] == "mol":
            # Check if molecule is a nitrile intermediate
            if checker.check_fg("Nitrile", node["smiles"]) and depth > 0:
                has_nitrile_intermediate = True
                print(f"Nitrile intermediate found at depth {depth}: {node['smiles']}")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Check if reaction converts nitrile to heterocycle
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has nitrile
            nitrile_reactants = [
                r for r in reactants if r and checker.check_fg("Nitrile", r)
            ]
            reactant_has_nitrile = len(nitrile_reactants) > 0

            # Comprehensive list of heterocyclic rings that can be formed from nitriles
            heterocycle_rings = [
                "pyrazole",
                "tetrazole",
                "triazole",
                "oxadiazole",
                "thiadiazole",
                "isoxazole",
                "isothiazole",
                "imidazole",
                "oxazole",
                "thiazole",
                "pyrimidine",
                "pyridazine",
                "pyrazine",
                "pyridine",
            ]

            # Check if product has any heterocyclic ring
            product_has_heterocycle = any(
                checker.check_ring(ring, product) for ring in heterocycle_rings
            )

            # Check if reactants already have the heterocycle
            reactants_combined = ".".join(reactants)
            reactants_have_heterocycle = any(
                checker.check_ring(ring, reactants_combined)
                for ring in heterocycle_rings
            )

            # Check for specific nitrile-to-heterocycle reactions
            nitrile_to_heterocycle_reactions = [
                "Azide-nitrile click cycloaddition to tetrazole",
                "Azide-nitrile click cycloaddition to triazole",
                "pyrazole",
                "tetrazole_terminal",
                "tetrazole_connect_regioisomere_1",
                "tetrazole_connect_regioisomere_2",
                "1,2,4-triazole_acetohydrazide",
                "1,2,4-triazole_carboxylic-acid/ester",
                "3-nitrile-pyridine",
                "oxadiazole",
                "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                "Huisgen 1,3 dipolar cycloaddition",
            ]

            reaction_is_nitrile_to_heterocycle = any(
                checker.check_reaction(rxn, rsmi)
                for rxn in nitrile_to_heterocycle_reactions
            )

            # If no specific reaction detected, check for general pattern of nitrile conversion to heterocycle
            if (
                reactant_has_nitrile
                and product_has_heterocycle
                and not reactants_have_heterocycle
            ):
                # This is a more general check for nitrile to heterocycle transformation
                if not reaction_is_nitrile_to_heterocycle:
                    # Check if the product has a new heterocycle not present in reactants
                    new_heterocycle = False
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product) and not checker.check_ring(
                            ring, reactants_combined
                        ):
                            new_heterocycle = True
                            print(f"New heterocycle {ring} detected in product")
                            break

                    reaction_is_nitrile_to_heterocycle = new_heterocycle

            # Check if this is a late-stage reaction (depth <= 3)
            if (
                reactant_has_nitrile
                and product_has_heterocycle
                and reaction_is_nitrile_to_heterocycle
                and depth <= 3
            ):
                nitrile_to_heterocycle = True
                print(f"Nitrile to heterocycle transformation found at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                print(
                    f"Nitrile reactant: {nitrile_reactants[0] if nitrile_reactants else 'None'}"
                )
                print(f"Heterocycle product: {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = has_nitrile_intermediate and nitrile_to_heterocycle
    print(f"Has nitrile intermediate: {has_nitrile_intermediate}")
    print(f"Has nitrile to heterocycle transformation: {nitrile_to_heterocycle}")
    print(f"Final result: {result}")

    return result
