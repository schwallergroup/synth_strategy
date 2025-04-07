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
    This function detects if the synthetic route involves late-stage heterocycle formation
    (cyclization in the final steps of the synthesis).
    """
    has_late_heterocycle = False

    # List of common heterocycles to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
    ]

    # List of common heterocycle formation reactions
    heterocycle_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Pictet-Spengler",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "1,2,4-triazole_acetohydrazide",
        "pyrazole",
        "oxadiazole",
        "Formation of NOS Heterocycles",
        "benzimidazole formation from aldehyde",
        "benzimidazole formation from acyl halide",
        "benzimidazole formation from ester/carboxylic acid",
        "benzimidazole formation (intramolecular)",
    ]

    # Additional functional groups that might indicate heterocycle formation
    heterocycle_related_fgs = [
        "Primary amine",
        "Secondary amine",
        "Aldehyde",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Azide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal has_late_heterocycle

        if node["type"] == "reaction" and depth <= 2:  # Check late-stage reactions (low depth)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    parts = rsmi.split(">")
                    reactants = parts[0]
                    product = parts[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check if this is a known heterocycle formation reaction
                    for rxn_type in heterocycle_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found heterocycle formation reaction: {rxn_type}")
                            has_late_heterocycle = True
                            return

                    # Check for heterocycle formation by comparing reactants and products
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and all(m for m in reactant_mols):
                        # Check if any heterocycle is present in the product but not in the reactants
                        for heterocycle in heterocycles:
                            product_has_heterocycle = checker.check_ring(heterocycle, product)
                            reactants_have_heterocycle = any(
                                checker.check_ring(heterocycle, r) for r in reactants.split(".")
                            )

                            if product_has_heterocycle and not reactants_have_heterocycle:
                                print(f"Found new heterocycle formation: {heterocycle}")
                                has_late_heterocycle = True
                                return

                        # Check for specific heterocycle patterns that might be missed
                        if (
                            "n1c2c(nc1)" in product
                            or "n1cnc2c1" in product
                            or "[nH]1c2c(nc1)" in product
                        ):
                            print("Found imidazole-type heterocycle formation pattern")
                            has_late_heterocycle = True
                            return

                        # Check for functional group changes that might indicate heterocycle formation
                        for fg in heterocycle_related_fgs:
                            reactants_have_fg = any(
                                checker.check_fg(fg, r) for r in reactants.split(".")
                            )
                            product_has_fg = checker.check_fg(fg, product)

                            if reactants_have_fg and not product_has_fg:
                                print(
                                    f"Functional group {fg} consumed in reaction - possible heterocycle formation"
                                )

                                # Additional check for new ring formation
                                try:
                                    reactant_rings = sum(
                                        len(m.GetRingInfo().AtomRings()) for m in reactant_mols if m
                                    )
                                    product_rings = len(product_mol.GetRingInfo().AtomRings())
                                    print(
                                        f"Ring structures in reactants: {reactant_rings}, in product: {product_rings}"
                                    )

                                    # Check if the product has any nitrogen-containing rings
                                    if "n" in product.lower() or "N" in product:
                                        print(
                                            "Product contains nitrogen in ring structure - likely heterocycle formation"
                                        )
                                        has_late_heterocycle = True
                                        return
                                except Exception as e:
                                    print(f"Error in ring analysis: {e}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late-stage heterocycle formation detected: {has_late_heterocycle}")
    return has_late_heterocycle
