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
    This function detects a strategy involving modification of a heterocyclic ring
    (e.g., dehydrogenation of piperidine, oxidation of heterocycles, ring opening/closing).
    """
    has_heterocycle_modification = False

    # Define heterocycles to check
    heterocycles = [
        "piperidine",
        "pyrrolidine",
        "morpholine",
        "thiomorpholine",
        "pyrrole",
        "furan",
        "thiophene",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "piperazine",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "aziridine",
        "azetidine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # Define reaction types that could modify heterocycles
    modification_reactions = [
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Hydrogenation (double to single)",
        "Hydrogenation (triple to double)",
        "Arene hydrogenation",
        "Aromatic hydroxylation",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_modification

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Method 1: Check for heterocycle transformation (one heterocycle to another)
                for reactant in reactants:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            # Check if product contains a different heterocycle
                            if any(
                                checker.check_ring(h, product) for h in heterocycles
                            ) and not checker.check_ring(heterocycle, product):
                                has_heterocycle_modification = True
                                print(
                                    f"Found heterocycle transformation from {heterocycle} at depth {depth}"
                                )
                                break

                # Method 2: Check for specific reaction types that modify heterocycles
                for reaction_type in modification_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        # Verify that a heterocycle is involved in the reaction
                        for reactant in reactants:
                            if any(
                                checker.check_ring(heterocycle, reactant)
                                for heterocycle in heterocycles
                            ):
                                has_heterocycle_modification = True
                                print(
                                    f"Found heterocycle modification via {reaction_type} at depth {depth}"
                                )
                                break

                # Method 3: Check for ring opening/closing reactions
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Check if number of rings changed
                        reactant_rings = reactant_mol.GetRingInfo().NumRings()
                        product_rings = product_mol.GetRingInfo().NumRings()

                        if reactant_rings != product_rings:
                            # Verify that a heterocycle is involved
                            if any(
                                checker.check_ring(heterocycle, reactant)
                                for heterocycle in heterocycles
                            ) or any(
                                checker.check_ring(heterocycle, product)
                                for heterocycle in heterocycles
                            ):
                                has_heterocycle_modification = True
                                print(f"Found heterocycle ring opening/closing at depth {depth}")

                # Method 4: Check for dehydrogenation of saturated heterocycles
                saturated_heterocycles = [
                    "piperidine",
                    "pyrrolidine",
                    "morpholine",
                    "thiomorpholine",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                ]
                unsaturated_heterocycles = ["pyridine", "pyrrole", "furan", "thiophene"]

                for i, reactant in enumerate(reactants):
                    for sat_heterocycle in saturated_heterocycles:
                        if checker.check_ring(sat_heterocycle, reactant):
                            for unsat_heterocycle in unsaturated_heterocycles:
                                if checker.check_ring(unsat_heterocycle, product):
                                    has_heterocycle_modification = True
                                    print(
                                        f"Found heterocycle dehydrogenation from {sat_heterocycle} to {unsat_heterocycle} at depth {depth}"
                                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if has_heterocycle_modification:
        print("Detected heterocycle modification strategy")

    return has_heterocycle_modification
