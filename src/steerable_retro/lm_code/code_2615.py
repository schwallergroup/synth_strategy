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
    Detects if the synthesis route involves functionalization of heterocycles
    (pyrazole, benzimidazole, and other nitrogen-containing heterocycles).
    """
    heterocycle_functionalization = False

    # List of heterocycles to check
    heterocycles = [
        "pyrazole",
        "benzimidazole",
        "imidazole",
        "triazole",
        "tetrazole",
        "oxazole",
        "thiazole",
        "isoxazole",
        "pyrrole",
    ]

    # List of functionalization reaction types
    functionalization_reactions = [
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Suzuki coupling with boronic acids",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Heck terminal vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Michael addition",
        "Methylation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_functionalization

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a functionalization reaction
            is_functionalization_rxn = False
            for rxn_type in functionalization_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    is_functionalization_rxn = True
                    print(f"Found functionalization reaction: {rxn_type}")
                    break

            if is_functionalization_rxn:
                # Check if product contains a heterocycle
                product_has_heterocycle = False
                product_heterocycle = None
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        product_has_heterocycle = True
                        product_heterocycle = heterocycle
                        print(f"Product contains {heterocycle}")
                        break

                if product_has_heterocycle:
                    # Check if any reactant also contains the same heterocycle
                    for reactant in reactants:
                        if checker.check_ring(product_heterocycle, reactant):
                            print(f"Reactant also contains {product_heterocycle}")
                            print(f"Heterocycle functionalization detected at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            heterocycle_functionalization = True
                            break

            # Alternative approach: check for heterocycle in both reactant and product
            # even if we don't recognize the specific reaction type
            if not heterocycle_functionalization:
                for heterocycle in heterocycles:
                    # Check if product contains the heterocycle
                    if checker.check_ring(heterocycle, product):
                        # Check if any reactant also contains the heterocycle
                        for reactant in reactants:
                            if checker.check_ring(heterocycle, reactant):
                                # Check if there's a change in functional groups
                                product_mol = Chem.MolFromSmiles(product)
                                reactant_mol = Chem.MolFromSmiles(reactant)

                                if product_mol and reactant_mol:
                                    # Count atoms as a simple way to detect modification
                                    if product_mol.GetNumAtoms() != reactant_mol.GetNumAtoms():
                                        print(
                                            f"Heterocycle {heterocycle} modified (atom count changed)"
                                        )
                                        print(f"Reaction SMILES: {rsmi}")
                                        heterocycle_functionalization = True
                                        break

                                    # Check for specific functional group additions
                                    for fg in [
                                        "Nitrile",
                                        "Ester",
                                        "Amide",
                                        "Halide",
                                        "Nitro group",
                                    ]:
                                        if (
                                            checker.check_fg(fg, product)
                                            and not checker.check_fg(fg, reactant)
                                        ) or (
                                            not checker.check_fg(fg, product)
                                            and checker.check_fg(fg, reactant)
                                        ):
                                            print(
                                                f"Heterocycle {heterocycle} functionalization detected with {fg} change"
                                            )
                                            print(f"Reaction SMILES: {rsmi}")
                                            heterocycle_functionalization = True
                                            break

                        if heterocycle_functionalization:
                            break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_functionalization
