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
    This function detects a linear synthesis strategy that builds heterocycles
    through sequential functional group transformations and ring formations.
    """
    # Track key transformations
    functional_group_mods = []
    heterocycle_formations = []
    reaction_count = 0

    # List of common heterocycles to check
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
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "thiophene",
        "triazole",
        "tetrazole",
        "oxadiazole",
        "thiadiazole",
        "isoxazole",
        "isothiazole",
        "benzotriazole",
        "tetrahydrofuran",
        "tetrahydropyran",
        "pyrrolidine",
        "piperidine",
    ]

    # List of functional group transformations to check
    fg_transformations = [
        ("Ether", "Phenol", "O-demethylation"),
        ("Nitro group", "Primary amine", "Nitro reduction"),
        ("Ester", "Carboxylic acid", "Ester hydrolysis"),
        ("Nitrile", "Primary amide", "Nitrile hydration"),
        ("Nitrile", "Primary amine", "Nitrile reduction"),
        ("Carboxylic acid", "Ester", "Esterification"),
        ("Carboxylic acid", "Primary amide", "Amidation"),
        ("Acyl halide", "Primary amide", "Amidation"),
        ("Aldehyde", "Primary alcohol", "Aldehyde reduction"),
        ("Ketone", "Secondary alcohol", "Ketone reduction"),
        ("Alkyne", "Alkene", "Alkyne reduction"),
        ("Primary halide", "Primary alcohol", "Halide hydrolysis"),
        ("Primary halide", "Primary amine", "Amination"),
        ("Azide", "Primary amine", "Azide reduction"),
    ]

    # List of ring-forming reactions
    ring_forming_reactions = [
        "Formation of NOS Heterocycles",
        "Paal-Knorr pyrrole synthesis",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "Huisgen_Ru-catalyzed_1,5_subst",
        "Huisgen_disubst-alkyne",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Azide-nitrile click cycloaddition to tetrazole",
        "Azide-nitrile click cycloaddition to triazole",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal functional_group_mods, heterocycle_formations, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for ring formation
                reactant_ring_counts = []
                reactant_has_heterocycle = False
                reactant_heterocycles = set()

                for reactant in reactants:
                    if not reactant:
                        continue

                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Count rings in reactant
                        ring_count = reactant_mol.GetRingInfo().NumRings()
                        reactant_ring_counts.append(ring_count)
                        print(f"  Reactant {reactant} has {ring_count} rings")

                        # Check if reactant already has a heterocycle
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                reactant_has_heterocycle = True
                                reactant_heterocycles.add(heterocycle)
                                print(f"  Reactant contains heterocycle: {heterocycle}")

                # Check product rings and heterocycles
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    product_ring_count = product_mol.GetRingInfo().NumRings()
                    product_heterocycles = set()

                    print(f"  Product {product} has {product_ring_count} rings")

                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product):
                            product_heterocycles.add(heterocycle)
                            print(f"  Product contains heterocycle: {heterocycle}")

                    # Detect heterocycle formation
                    new_heterocycles = product_heterocycles - reactant_heterocycles
                    if new_heterocycles:
                        print(
                            f"  New heterocycle(s) formed: {', '.join(new_heterocycles)}"
                        )
                        heterocycle_formations.append((depth, list(new_heterocycles)))

                    # Detect ring formation (even if not heterocycle)
                    max_reactant_rings = (
                        max(reactant_ring_counts) if reactant_ring_counts else 0
                    )
                    if product_ring_count > max_reactant_rings:
                        print(
                            f"  Ring formation detected: {max_reactant_rings} → {product_ring_count}"
                        )
                        heterocycle_formations.append((depth, ["ring_formation"]))

                # Check for functional group modifications
                for start_fg, end_fg, name in fg_transformations:
                    reactant_has_start_fg = any(
                        checker.check_fg(start_fg, r) for r in reactants if r
                    )
                    product_has_end_fg = checker.check_fg(end_fg, product)

                    if reactant_has_start_fg and product_has_end_fg:
                        print(
                            f"  Functional group transformation detected: {name} ({start_fg} → {end_fg})"
                        )
                        functional_group_mods.append((name, depth))

                # Specifically check for O-demethylation (methoxy to hydroxyl conversion)
                for reactant in reactants:
                    if not reactant:
                        continue
                    if checker.check_fg("Ether", reactant) and checker.check_fg(
                        "Phenol", product
                    ):
                        print(
                            f"  Functional group transformation detected: O-demethylation (Methoxy → Hydroxyl)"
                        )
                        functional_group_mods.append(("O-demethylation", depth))
                        break

                # Check for specific reactions
                for reaction_type in ring_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"  Reaction detected: {reaction_type}")
                        heterocycle_formations.append((depth, [reaction_type]))
                        functional_group_mods.append((reaction_type, depth))

                # Check for other common reactions that modify functional groups
                common_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Reduction of nitro groups to amines",
                    "Reduction of nitrile to amine",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Reduction of ester to primary alcohol",
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of nitrile to carboxylic acid",
                    "Oxidation of amide to carboxylic acid",
                    "Cleavage of methoxy ethers to alcohols",
                    "Ether cleavage to primary alcohol",
                ]

                for reaction_type in common_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"  Reaction detected: {reaction_type}")
                        functional_group_mods.append((reaction_type, depth))

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Remove duplicates while preserving order
    unique_fg_mods = []
    seen_fg_mods = set()
    for mod, depth in functional_group_mods:
        if mod not in seen_fg_mods:
            unique_fg_mods.append((mod, depth))
            seen_fg_mods.add(mod)

    unique_heterocycle_formations = []
    seen_heterocycle_formations = set()
    for depth, formations in heterocycle_formations:
        for formation in formations:
            if formation not in seen_heterocycle_formations:
                unique_heterocycle_formations.append((formation, depth))
                seen_heterocycle_formations.add(formation)

    # Check if we have a linear synthesis with heterocycle construction
    has_functional_group_mods = len(unique_fg_mods) >= 1
    has_heterocycle_formation = len(unique_heterocycle_formations) >= 1
    is_linear = reaction_count >= 3  # At least 3 reactions in sequence

    print(
        f"FG mods: {len(unique_fg_mods)}, Heterocycle formations: {len(unique_heterocycle_formations)}, Reactions: {reaction_count}"
    )
    print(f"Functional group modifications: {unique_fg_mods}")
    print(f"Heterocycle formations: {unique_heterocycle_formations}")

    return has_functional_group_mods and has_heterocycle_formation and is_linear
