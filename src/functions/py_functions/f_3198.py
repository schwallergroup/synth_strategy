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
    This function detects a linear synthesis strategy with sequential
    functionalization of a heterocycle scaffold.
    """
    # Track key features
    reaction_count = 0
    linear_steps = 0
    has_heterocycle_core = False
    functional_group_changes = 0

    # List of heterocyclic rings to check
    heterocycle_rings = [
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
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "thiophene",
    ]

    # List of functional groups to track
    functional_groups = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Aldehyde",
        "Ketone",
        "Phenol",
        "Nitro group",
        "Boronic acid",
        "Boronic ester",
    ]

    # Track functional groups present in each molecule
    molecule_functional_groups = {}

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_steps, has_heterocycle_core, functional_group_changes

        if node["type"] == "mol":
            if node["smiles"]:
                mol_smiles = node["smiles"]
                molecule_functional_groups[mol_smiles] = set()

                # Check for heterocycle core in molecules
                heterocycle_found = False
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, mol_smiles):
                        has_heterocycle_core = True
                        heterocycle_found = True
                        print(
                            f"Heterocycle {ring} detected in molecule at depth {depth}"
                        )
                        break

                # If heterocycle found, check for functional groups
                if heterocycle_found:
                    for fg in functional_groups:
                        if checker.check_fg(fg, mol_smiles):
                            molecule_functional_groups[mol_smiles].add(fg)
                            print(
                                f"Functional group {fg} detected in molecule at depth {depth}"
                            )

        elif node["type"] == "reaction":
            reaction_count += 1

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a linear step (1-2 reactants)
                if len(reactants) <= 2:
                    linear_steps += 1
                    print(
                        f"Linear step detected at depth {depth} with {len(reactants)} reactants"
                    )

                # Check for heterocycle in reactants and products
                product_has_heterocycle = False
                reactants_with_heterocycle = []

                for ring in heterocycle_rings:
                    if checker.check_ring(ring, product):
                        has_heterocycle_core = True
                        product_has_heterocycle = True
                        print(
                            f"Heterocycle {ring} detected in product at depth {depth}"
                        )
                        break

                for reactant in reactants:
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, reactant):
                            has_heterocycle_core = True
                            reactants_with_heterocycle.append(reactant)
                            print(
                                f"Heterocycle {ring} detected in reactant at depth {depth}"
                            )
                            break

                # Check for functional group changes
                # First check predefined reaction types
                functionalization_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation of alcohol to carboxylic acid",
                    "Reduction of ester to primary alcohol",
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Alkylation of amines",
                    "Suzuki coupling with boronic acids",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Suzuki coupling with boronic esters",
                    "Heck terminal vinyl",
                    "Sonogashira alkyne_aryl halide",
                    "Williamson Ether Synthesis",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                ]

                reaction_detected = False
                for reaction_type in functionalization_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        reaction_detected = True
                        functional_group_changes += 1
                        print(
                            f"Functional group change detected: {reaction_type} at depth {depth}"
                        )
                        break

                # If no predefined reaction detected, check for functional group changes directly
                if (
                    not reaction_detected
                    and product_has_heterocycle
                    and reactants_with_heterocycle
                ):
                    # Check functional groups in product
                    product_fgs = set()
                    for fg in functional_groups:
                        if checker.check_fg(fg, product):
                            product_fgs.add(fg)

                    # Check functional groups in heterocycle-containing reactants
                    for reactant in reactants_with_heterocycle:
                        reactant_fgs = set()
                        for fg in functional_groups:
                            if checker.check_fg(fg, reactant):
                                reactant_fgs.add(fg)

                        # If there's a difference in functional groups, count it as a change
                        if product_fgs != reactant_fgs:
                            functional_group_changes += 1
                            print(
                                f"Functional group change detected between reactant and product at depth {depth}"
                            )
                            print(f"  Reactant FGs: {reactant_fgs}")
                            print(f"  Product FGs: {product_fgs}")
                            break

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this strategy is present
    strategy_present = (
        has_heterocycle_core
        and reaction_count >= 3
        and linear_steps >= 3
        and functional_group_changes >= 1
    )

    print(f"Linear heterocycle functionalization strategy detection results:")
    print(f"  Heterocycle core: {has_heterocycle_core}")
    print(f"  Total reactions: {reaction_count}")
    print(f"  Linear steps: {linear_steps}")
    print(f"  Functional group changes: {functional_group_changes}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
