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
    This function detects the construction of nitrogen-rich heterocycles in the synthesis route.
    """
    n_heterocycle_construction_found = False

    def dfs_traverse(node):
        nonlocal n_heterocycle_construction_found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # List of nitrogen-rich heterocycles to check for
                n_rich_rings = [
                    "triazole",
                    "tetrazole",
                    "imidazole",
                    "pyrazole",
                    "oxadiazole",
                    "thiadiazole",
                    "pyrimidine",
                    "pyrazine",
                    "pyridazine",
                    "benzotriazole",
                    "indazole",
                    "benzimidazole",
                ]

                # Check if product contains any nitrogen-rich heterocycles
                product_has_n_rich_ring = False
                n_rich_ring_found = None

                for ring in n_rich_rings:
                    if checker.check_ring(ring, product):
                        product_has_n_rich_ring = True
                        n_rich_ring_found = ring
                        print(f"Found {ring} in product: {product}")
                        break

                if product_has_n_rich_ring and n_rich_ring_found:
                    # Check if any reactant contains the same nitrogen-rich heterocycle
                    reactants_have_same_ring = False
                    for reactant in reactants:
                        if checker.check_ring(n_rich_ring_found, reactant):
                            reactants_have_same_ring = True
                            print(f"Found {n_rich_ring_found} in reactant: {reactant}")
                            break

                    if not reactants_have_same_ring:
                        # Check if this is a known reaction type for forming nitrogen-rich heterocycles
                        n_heterocycle_rxn_types = [
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                            "Huisgen 1,3 dipolar cycloaddition",
                            "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                            "Huisgen_Cu-catalyzed_1,4-subst",
                            "Huisgen_Ru-catalyzed_1,5_subst",
                            "Huisgen_disubst-alkyne",
                            "Azide-nitrile click cycloaddition to tetrazole",
                            "Azide-nitrile click cycloaddition to triazole",
                            "tetrazole_terminal",
                            "tetrazole_connect_regioisomere_1",
                            "tetrazole_connect_regioisomere_2",
                            "1,2,4-triazole_acetohydrazide",
                            "1,2,4-triazole_carboxylic-acid/ester",
                            "Pyrazole formation",
                            "pyrazole",
                            "imidazole",
                            "triaryl-imidazole",
                            "oxadiazole",
                            "benzimidazole_derivatives_carboxylic-acid/ester",
                            "benzimidazole_derivatives_aldehyde",
                            "benzotriazole",
                        ]

                        rxn_match_found = False
                        for rxn_type in n_heterocycle_rxn_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                n_heterocycle_construction_found = True
                                rxn_match_found = True
                                print(
                                    f"Nitrogen-rich heterocycle construction detected via reaction: {rxn_type}"
                                )
                                break

                        # If no specific reaction type was found but we still have a new nitrogen-rich heterocycle,
                        # it's still a construction
                        if not rxn_match_found:
                            n_heterocycle_construction_found = True
                            print(
                                f"Nitrogen-rich heterocycle construction detected: {n_rich_ring_found}"
                            )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_heterocycle_construction_found
