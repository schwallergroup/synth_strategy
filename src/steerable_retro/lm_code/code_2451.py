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
    This function detects if the synthesis involves modifications to a heterocyclic aromatic core
    like a benzoxazole-like structure.
    """
    # Track if we've found a heterocyclic core modification
    heterocyclic_cores_modified = False

    # List of heterocyclic aromatic cores to check
    heterocyclic_rings = [
        "benzoxazole",
        "benzimidazole",
        "benzothiazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "furan",
        "pyrrole",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
    ]

    # Reactions that commonly modify heterocyclic cores
    core_modifying_reactions = [
        "Friedel-Crafts acylation",
        "Friedel-Crafts alkylation",
        "Aromatic bromination",
        "Aromatic chlorination",
        "Aromatic fluorination",
        "Aromatic iodination",
        "Aromatic nitration",
        "N-arylation",
        "Buchwald-Hartwig",
        "Suzuki coupling",
        "Heck reaction",
        "Sonogashira",
        "Minisci",
        "Directed ortho metalation of arenes",
        "Acylation of Nitrogen Nucleophiles",
        "Esterification",
        "Williamson Ether Synthesis",
        "O-alkylation",
        "S-alkylation",
        "N-alkylation",
        "Oxidation",
        "Reduction",
    ]

    # Functional groups that might be added or modified on heterocyclic cores
    modifying_functional_groups = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Phenol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aniline",
        "Carboxylic acid",
        "Ester",
        "Amide",
        "Nitrile",
        "Aldehyde",
        "Ketone",
        "Halide",
        "Nitro group",
        "Ether",
        "Thiol",
        "Sulfide",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocyclic_cores_modified

        # Check reaction nodes for core modifications
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains a heterocyclic core
            reactant_cores = {}
            for reactant in reactants:
                for ring in heterocyclic_rings:
                    if checker.check_ring(ring, reactant):
                        if ring not in reactant_cores:
                            reactant_cores[ring] = []
                        reactant_cores[ring].append(reactant)
                        print(f"Found {ring} in reactant: {reactant}")

            # Check if product contains a heterocyclic core
            product_cores = []
            for ring in heterocyclic_rings:
                if checker.check_ring(ring, product):
                    product_cores.append(ring)
                    print(f"Found {ring} in product: {product}")

            # If both reactants and product have heterocyclic cores
            if reactant_cores and product_cores:
                # First check if this is a known core-modifying reaction
                for reaction_type in core_modifying_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected core-modifying reaction: {reaction_type}")

                        # Now verify that functional groups on the core are changing
                        # by checking for functional groups in reactants and product
                        reactant_fgs = set()
                        product_fgs = set()

                        for fg in modifying_functional_groups:
                            # Check reactants for functional groups
                            for reactant in reactants:
                                if checker.check_fg(fg, reactant):
                                    reactant_fgs.add(fg)
                                    print(f"Found {fg} in reactant")

                            # Check product for functional groups
                            if checker.check_fg(fg, product):
                                product_fgs.add(fg)
                                print(f"Found {fg} in product")

                        # If there's a difference in functional groups, it's a core modification
                        if reactant_fgs != product_fgs:
                            heterocyclic_cores_modified = True
                            print(
                                f"Confirmed heterocyclic core modification: functional groups changed from {reactant_fgs} to {product_fgs}"
                            )
                            break

                # Check for specific functional group transformations on the core
                if not heterocyclic_cores_modified:
                    # Check for common transformations like hydroxylation, alkylation, etc.
                    transformations = [
                        ("Phenol", "Ether"),  # Hydroxyl to ether
                        ("Primary alcohol", "Ester"),  # Alcohol to ester
                        ("Carboxylic acid", "Ester"),  # Acid to ester
                        ("Carboxylic acid", "Amide"),  # Acid to amide
                        ("Aldehyde", "Primary alcohol"),  # Aldehyde reduction
                        ("Ketone", "Secondary alcohol"),  # Ketone reduction
                        ("Primary amine", "Secondary amine"),  # Amine alkylation
                        ("Secondary amine", "Tertiary amine"),  # Amine alkylation
                    ]

                    for from_fg, to_fg in transformations:
                        # Check if transformation occurs from reactant to product
                        reactant_has_from = any(checker.check_fg(from_fg, r) for r in reactants)
                        product_has_to = checker.check_fg(to_fg, product)

                        if reactant_has_from and product_has_to:
                            print(f"Detected transformation from {from_fg} to {to_fg}")
                            heterocyclic_cores_modified = True
                            break

                # If we still haven't confirmed a modification, check if the core structure is preserved
                # but with different substituents (using atom mapping)
                if not heterocyclic_cores_modified:
                    # Check if the reaction involves a modification to the heterocyclic core
                    # by looking at specific functional group changes on the core
                    for ring in set(reactant_cores.keys()).intersection(set(product_cores)):
                        # If the same ring type is in both reactants and product, check for modifications
                        print(f"Same ring type {ring} found in both reactants and product")

                        # Check for specific functional group changes on the heterocyclic core
                        for fg in modifying_functional_groups:
                            reactant_has_fg = any(checker.check_fg(fg, r) for r in reactants)
                            product_has_fg = checker.check_fg(fg, product)

                            if reactant_has_fg != product_has_fg:
                                print(f"Functional group {fg} changed on the heterocyclic core")
                                heterocyclic_cores_modified = True
                                break

                # Last resort: check if any reaction is occurring on the heterocyclic core
                # by looking at the reaction type
                if not heterocyclic_cores_modified:
                    # Check for any reaction type that might modify a heterocyclic core
                    for reaction_type in [
                        "Acylation",
                        "Alkylation",
                        "Amination",
                        "Esterification",
                        "Etherification",
                        "Halogenation",
                        "Hydroxylation",
                        "Nitration",
                        "Oxidation",
                        "Reduction",
                        "Sulfonation",
                    ]:
                        # Use a partial match approach since reaction_type might be part of a longer name
                        for full_reaction in core_modifying_reactions:
                            if reaction_type in full_reaction and checker.check_reaction(
                                full_reaction, rsmi
                            ):
                                print(
                                    f"Detected potential core-modifying reaction: {full_reaction}"
                                )
                                heterocyclic_cores_modified = True
                                break

                        if heterocyclic_cores_modified:
                            break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Synthesis {'involves' if heterocyclic_cores_modified else 'does not involve'} heterocyclic aromatic core modification"
    )
    return heterocyclic_cores_modified
