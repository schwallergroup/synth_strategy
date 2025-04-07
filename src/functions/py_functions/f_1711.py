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
    Detects a strategy involving late-stage heterocycle modification,
    specifically the conversion of a saturated heterocycle to an aromatic one.
    In retrosynthetic analysis, this means finding an aromatic heterocycle
    that is converted to a saturated one.
    """
    found_late_stage_heterocycle_mod = False

    # Define pairs of saturated and aromatic heterocycles to check
    heterocycle_pairs = [
        # Nitrogen-containing
        ("piperidine", "pyridine"),
        ("pyrrolidine", "pyrrole"),
        ("piperazine", "pyrazine"),
        # Oxygen-containing
        ("tetrahydrofuran", "furan"),
        ("tetrahydropyran", "pyran"),
        ("dioxane", "dioxene"),
        # Sulfur-containing
        ("thiomorpholine", "thiopyran"),
        ("thiolane", "thiophene"),
        ("thiane", "thiopyran"),
        # Mixed heterocycles
        ("oxazolidine", "oxazole"),
        ("thiazolidine", "thiazole"),
        ("imidazolidine", "imidazole"),
        ("pyrroline", "pyrrole"),
        ("dihydropyridine", "pyridine"),
        ("dihydropyrazine", "pyrazine"),
        ("dihydropyrimidine", "pyrimidine"),
    ]

    # List of oxidation/dehydrogenation reaction types
    oxidation_reactions = [
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Arene hydrogenation",  # This works in reverse for dehydrogenation
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation of amide to carboxylic acid",
        "Oxidation of alkene to aldehyde",
        "Oxidative esterification of primary alcohols",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_heterocycle_mod

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check reactions up to depth 2 (late stage)
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, the product of the reaction is actually the reactant in forward synthesis
                # and the reactants are the products in forward synthesis
                for reactant in reactants:
                    for sat_ring, arom_ring in heterocycle_pairs:
                        # Check if reactant (product in forward synthesis) contains aromatic ring
                        # and product (reactant in forward synthesis) contains saturated ring
                        if checker.check_ring(
                            arom_ring, reactant
                        ) and checker.check_ring(sat_ring, product):
                            print(
                                f"Found potential heterocycle pair: {arom_ring} in reactant, {sat_ring} in product"
                            )

                            # Check for oxidation/dehydrogenation reaction types
                            # In retrosynthesis, we're looking for the reverse of these reactions
                            for rxn_type in oxidation_reactions:
                                if checker.check_reaction(rxn_type, rsmi):
                                    print(
                                        f"Found late-stage heterocycle modification ({arom_ring} to {sat_ring}) with reaction type: {rxn_type}"
                                    )
                                    found_late_stage_heterocycle_mod = True
                                    # Don't return, continue checking for other instances

                            # If no specific reaction type matched but we have the right structural change
                            if not found_late_stage_heterocycle_mod:
                                # Check for general oxidation/dehydrogenation patterns
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                product_mol = Chem.MolFromSmiles(product)

                                if reactant_mol and product_mol:
                                    # Count hydrogens in reactant vs product
                                    reactant_h_count = sum(
                                        [
                                            atom.GetTotalNumHs()
                                            for atom in reactant_mol.GetAtoms()
                                        ]
                                    )
                                    product_h_count = sum(
                                        [
                                            atom.GetTotalNumHs()
                                            for atom in product_mol.GetAtoms()
                                        ]
                                    )

                                    # In retrosynthesis, the product (reactant in forward direction)
                                    # should have more hydrogens than the reactant (product in forward direction)
                                    if (
                                        reactant_h_count < product_h_count
                                        and depth <= 1
                                    ):
                                        print(
                                            f"Found late-stage heterocycle modification ({arom_ring} to {sat_ring}) based on hydrogen count difference"
                                        )
                                        print(
                                            f"Hydrogen count: reactant={reactant_h_count}, product={product_h_count}"
                                        )
                                        found_late_stage_heterocycle_mod = True
                                        # Don't return, continue checking for other instances

                                    # Check specifically for arene hydrogenation in reverse
                                    if (
                                        checker.check_reaction(
                                            "Arene hydrogenation", rsmi
                                        )
                                        and depth <= 1
                                    ):
                                        # In retrosynthesis, this would be dehydrogenation in forward direction
                                        print(
                                            f"Found late-stage heterocycle modification ({arom_ring} to {sat_ring}) via arene hydrogenation (reverse)"
                                        )
                                        found_late_stage_heterocycle_mod = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage heterocycle modification strategy present: {found_late_stage_heterocycle_mod}"
    )
    return found_late_stage_heterocycle_mod
