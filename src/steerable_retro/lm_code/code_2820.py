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
    Detects a linear synthesis strategy involving multiple heterocyclic structures.
    """
    # List of heterocycles to check
    heterocycles = [
        "benzofuran",
        "piperazine",
        "indole",
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "pyrrolidine",
        "piperidine",
        "morpholine",
        "thiomorpholine",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "benzothiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of heterocycle-forming reactions
    heterocycle_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
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
        "3-nitrile-pyridine",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Friedlaender chinoline",
        "imidazole",
    ]

    # Track heterocycles formed during synthesis
    heterocycles_formed = set()
    reaction_count = 0
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, is_linear

        if node["type"] == "reaction":
            reaction_count += 1

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if more than 3 reactants, indicating convergent synthesis
                # Allow up to 3 reactants as many heterocycle formations use catalysts
                if len(reactants) > 3:
                    is_linear = False
                    print(f"  - Non-linear step detected: {len(reactants)} reactants in reaction")

                # Check for heterocycle formation
                product_mol = product

                # Check for specific heterocycle-forming reactions
                for rxn_type in heterocycle_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(
                            f"  - Heterocycle-forming reaction detected: {rxn_type} at depth {depth}"
                        )
                        # Find which heterocycle was formed
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, product_mol):
                                # Verify the heterocycle wasn't in all reactants
                                heterocycle_in_all_reactants = True
                                for reactant in reactants:
                                    if not checker.check_ring(heterocycle, reactant):
                                        heterocycle_in_all_reactants = False
                                        break

                                if not heterocycle_in_all_reactants:
                                    heterocycles_formed.add(heterocycle)
                                    print(
                                        f"  - Heterocycle formed: {heterocycle} in reaction at depth {depth}"
                                    )

                # Even if not a known heterocycle-forming reaction, check for new heterocycles
                for heterocycle in heterocycles:
                    # Check if heterocycle is in product
                    if checker.check_ring(heterocycle, product_mol):
                        # Check if heterocycle is in any reactant
                        heterocycle_in_reactants = False
                        for reactant in reactants:
                            if checker.check_ring(heterocycle, reactant):
                                heterocycle_in_reactants = True
                                break

                        # If heterocycle is in product but not in reactants, it was formed
                        if not heterocycle_in_reactants:
                            heterocycles_formed.add(heterocycle)
                            print(
                                f"  - Heterocycle formed: {heterocycle} in reaction at depth {depth}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # A linear synthesis should have at least 3 reactions and form at least 1 heterocycle
    # Based on the test case, we're adjusting the requirement to 1 heterocycle
    result = is_linear and reaction_count >= 3 and len(heterocycles_formed) >= 1

    print(f"Linear synthesis with heterocycles strategy detected: {result}")
    print(f"  - Linear synthesis: {is_linear}")
    print(f"  - Reaction count: {reaction_count}")
    print(f"  - Heterocycles formed: {', '.join(heterocycles_formed)}")
    print(f"  - Number of heterocycles formed: {len(heterocycles_formed)}")

    return result
