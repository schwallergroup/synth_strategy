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
    Detects if the synthesis route involves functionalization of heterocycles.
    """
    result = False

    # List of heterocycles to check
    heterocycles = [
        "furan",
        "pyran",
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
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "thiophene",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indazole",
        "benzotriazole",
    ]

    # Reactions commonly used for heterocycle functionalization
    functionalization_reactions = [
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Suzuki coupling",
        "Buchwald-Hartwig",
        "N-arylation",
        "Heck terminal vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Negishi coupling",
        "Directed ortho metalation of arenes",
        "Minisci (para)",
        "Minisci (ortho)",
        "Catellani reaction ortho",
        "Catellani reaction para",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycles in reactants and product
                heterocycle_in_reactants = False
                heterocycle_reactant = None
                heterocycle_found = None

                for reactant in reactants_smiles:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_in_reactants = True
                            heterocycle_reactant = reactant
                            heterocycle_found = heterocycle
                            print(f"Found {heterocycle} in reactant: {reactant}")
                            break
                    if heterocycle_in_reactants:
                        break

                heterocycle_in_product = False
                if heterocycle_found and checker.check_ring(heterocycle_found, product_smiles):
                    heterocycle_in_product = True
                    print(f"Found {heterocycle_found} in product: {product_smiles}")

                # Check if this is a known functionalization reaction
                is_functionalization_rxn = False
                for rxn_type in functionalization_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_functionalization_rxn = True
                        print(f"Detected heterocycle functionalization reaction: {rxn_type}")
                        break

                if heterocycle_in_reactants and heterocycle_in_product:
                    # The heterocycle is preserved in the reaction
                    if is_functionalization_rxn:
                        # It's a known functionalization reaction
                        print("Confirmed heterocycle functionalization strategy via known reaction")
                        result = True
                        return
                    else:
                        # Check if the heterocycle is modified by comparing molecules
                        reactant_mol = Chem.MolFromSmiles(heterocycle_reactant)
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if reactant_mol and product_mol:
                            # If molecules are different but contain the same heterocycle,
                            # it's likely a functionalization
                            if Chem.MolToSmiles(reactant_mol) != Chem.MolToSmiles(product_mol):
                                # Additional check: look for common functionalization patterns
                                if (
                                    (
                                        checker.check_fg("Aromatic halide", heterocycle_reactant)
                                        and not checker.check_fg("Aromatic halide", product_smiles)
                                    )
                                    or (
                                        not checker.check_fg("Nitrile", heterocycle_reactant)
                                        and checker.check_fg("Nitrile", product_smiles)
                                    )
                                    or (
                                        not checker.check_fg("Ester", heterocycle_reactant)
                                        and checker.check_fg("Ester", product_smiles)
                                    )
                                    or (
                                        not checker.check_fg("Ketone", heterocycle_reactant)
                                        and checker.check_fg("Ketone", product_smiles)
                                    )
                                    or (
                                        not checker.check_fg(
                                            "Carboxylic acid", heterocycle_reactant
                                        )
                                        and checker.check_fg("Carboxylic acid", product_smiles)
                                    )
                                ):
                                    print(
                                        f"Detected heterocycle functionalization via functional group change"
                                    )
                                    print(f"Reactant: {heterocycle_reactant}")
                                    print(f"Product: {product_smiles}")
                                    result = True
                                    return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
