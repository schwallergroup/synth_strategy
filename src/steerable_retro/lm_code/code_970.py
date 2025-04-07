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
    This function detects the synthesis of a nitrogen-rich molecule containing
    multiple nitrogen-containing heterocycles (pyridine, piperidine, benzimidazole, etc.).
    """
    # List of nitrogen-containing heterocycles to check
    n_heterocycles = [
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
        "piperazine",
        "morpholine",
        "aziridine",
        "azetidine",
        "azepane",
        "diazepane",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "acridine",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indazole",
        "benzotriazole",
    ]

    # Track if we found synthesis of nitrogen-rich heterocycles
    found_n_rich_synthesis = False

    # Track heterocycles in final product
    final_product_heterocycles = set()
    final_product_n_count = 0
    has_nitrile = False

    # Track heterocycles created in reactions
    synthesized_heterocycles = set()
    heterocycle_forming_reactions = False

    # Track heterocycles in starting materials
    starting_material_heterocycles = set()

    # Track heterocycle modifications
    modified_heterocycles = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_rich_synthesis, final_product_heterocycles, final_product_n_count
        nonlocal synthesized_heterocycles, has_nitrile, heterocycle_forming_reactions
        nonlocal starting_material_heterocycles, modified_heterocycles

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # For the final product (depth 0), check all nitrogen heterocycles
            if depth == 0:
                # Count nitrogen atoms
                try:
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetAtomicNum() == 7:  # Nitrogen
                                final_product_n_count += 1
                except Exception as e:
                    print(f"Error counting nitrogen atoms: {e}")

                # Check for nitrogen heterocycles
                for ring in n_heterocycles:
                    try:
                        if checker.check_ring(ring, mol_smiles):
                            final_product_heterocycles.add(ring)
                            print(f"Found {ring} in final product")
                    except Exception as e:
                        print(f"Error checking for {ring}: {e}")

                # Also check for cyano group (separate from heterocycles)
                try:
                    if checker.check_fg("Nitrile", mol_smiles):
                        has_nitrile = True
                        print(f"Found Nitrile in final product")
                except Exception as e:
                    print(f"Error checking for Nitrile: {e}")

            # Check if this is a starting material (in_stock)
            elif node.get("in_stock", False):
                # Check for nitrogen heterocycles in starting materials
                for ring in n_heterocycles:
                    try:
                        if checker.check_ring(ring, mol_smiles):
                            starting_material_heterocycles.add(ring)
                            print(f"Found {ring} in starting material")
                    except Exception as e:
                        print(f"Error checking for {ring} in starting material: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction creates any nitrogen heterocycles
            # by comparing reactants and products
            product_rings = set()
            reactant_rings = set()

            # Check rings in product
            for ring in n_heterocycles:
                try:
                    if checker.check_ring(ring, product):
                        product_rings.add(ring)
                except Exception as e:
                    print(f"Error checking for {ring} in product: {e}")

            # Check rings in reactants
            for reactant in reactants:
                for ring in n_heterocycles:
                    try:
                        if checker.check_ring(ring, reactant):
                            reactant_rings.add(ring)
                    except Exception as e:
                        print(f"Error checking for {ring} in reactant: {e}")

            # Rings in product but not in reactants were synthesized
            newly_synthesized = product_rings - reactant_rings
            for ring in newly_synthesized:
                synthesized_heterocycles.add(ring)
                print(f"Synthesized {ring} in reaction at depth {depth}")

            # Check if existing heterocycles were modified (different count or different types)
            if len(product_rings) != len(reactant_rings) or product_rings != reactant_rings:
                if len(reactant_rings) > 0 and len(product_rings) > 0:
                    modified_heterocycles = True
                    print(f"Modified heterocycles in reaction at depth {depth}")

            # Check for specific heterocycle-forming reactions
            try:
                # Expanded list of heterocycle-forming reactions
                if (
                    checker.check_reaction("benzimidazole_derivatives_aldehyde", rsmi)
                    or checker.check_reaction(
                        "benzimidazole_derivatives_carboxylic-acid/ester", rsmi
                    )
                    or checker.check_reaction("benzoxazole_arom-aldehyde", rsmi)
                    or checker.check_reaction("benzoxazole_carboxylic-acid", rsmi)
                    or checker.check_reaction("thiazole", rsmi)
                    or checker.check_reaction("tetrazole_terminal", rsmi)
                    or checker.check_reaction("tetrazole_connect_regioisomere_1", rsmi)
                    or checker.check_reaction("tetrazole_connect_regioisomere_2", rsmi)
                    or checker.check_reaction("pyrazole", rsmi)
                    or checker.check_reaction("Fischer indole", rsmi)
                    or checker.check_reaction("indole", rsmi)
                    or checker.check_reaction("imidazole", rsmi)
                    or checker.check_reaction("1,2,4-triazole_acetohydrazide", rsmi)
                    or checker.check_reaction("1,2,4-triazole_carboxylic-acid/ester", rsmi)
                    or checker.check_reaction("Paal-Knorr pyrrole", rsmi)
                    or checker.check_reaction("triaryl-imidazole", rsmi)
                    or checker.check_reaction("oxadiazole", rsmi)
                    or checker.check_reaction("Huisgen_Cu-catalyzed_1,4-subst", rsmi)
                    or checker.check_reaction("Huisgen_Ru-catalyzed_1,5_subst", rsmi)
                    or checker.check_reaction("Huisgen_disubst-alkyne", rsmi)
                    or checker.check_reaction("Pictet-Spengler", rsmi)
                    or checker.check_reaction("Niementowski_quinazoline", rsmi)
                    or checker.check_reaction("N-arylation", rsmi)
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Azide-nitrile click cycloaddition to tetrazole", rsmi
                    )
                    or checker.check_reaction("Azide-nitrile click cycloaddition to triazole", rsmi)
                    or checker.check_reaction(
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                ):

                    print(f"Detected heterocycle-forming reaction at depth {depth}")
                    heterocycle_forming_reactions = True
            except Exception as e:
                print(f"Error checking for heterocycle-forming reactions: {e}")

            # Check for reactions that modify nitrogen-containing rings
            try:
                if (
                    checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Goldberg coupling aryl amine-aryl chloride", rsmi)
                    or checker.check_reaction("Goldberg coupling aryl amide-aryl chloride", rsmi)
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                ):

                    # Check if any reactant has a nitrogen heterocycle
                    for reactant in reactants:
                        for ring in n_heterocycles:
                            if checker.check_ring(ring, reactant):
                                print(f"Detected heterocycle-modifying reaction at depth {depth}")
                                modified_heterocycles = True
                                break
            except Exception as e:
                print(f"Error checking for heterocycle-modifying reactions: {e}")

            # Check for reactions involving nitrile groups that might lead to heterocycles
            try:
                if (
                    checker.check_reaction("Nitrile to amide", rsmi)
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                ):

                    # Check if product has a nitrogen heterocycle
                    for ring in n_heterocycles:
                        if checker.check_ring(ring, product):
                            print(
                                f"Detected nitrile-involving heterocycle formation at depth {depth}"
                            )
                            heterocycle_forming_reactions = True
                            break
            except Exception as e:
                print(f"Error checking for nitrile-involving reactions: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if we found a nitrogen-rich heterocycle synthesis
    # Criteria:
    # 1. Final product has at least 3 different nitrogen heterocycles or
    #    at least 2 heterocycles and a cyano group
    # 2. Final product has at least 4 nitrogen atoms
    # 3. At least one of:
    #    a) A heterocycle was synthesized during the route
    #    b) A heterocycle-forming reaction was detected
    #    c) Existing heterocycles were modified
    #    d) Final product heterocycles are different from starting material heterocycles

    heterocycle_count = len(final_product_heterocycles)

    # Check if final product heterocycles are different from starting materials
    heterocycle_difference = final_product_heterocycles - starting_material_heterocycles
    heterocycle_transformation = len(heterocycle_difference) > 0

    print(f"Final evaluation:")
    print(f"- Heterocycle count: {heterocycle_count}")
    print(f"- Has nitrile: {has_nitrile}")
    print(f"- Nitrogen count: {final_product_n_count}")
    print(f"- Synthesized heterocycles: {synthesized_heterocycles}")
    print(f"- Heterocycle-forming reactions detected: {heterocycle_forming_reactions}")
    print(f"- Modified heterocycles: {modified_heterocycles}")
    print(f"- Starting material heterocycles: {starting_material_heterocycles}")
    print(f"- Final product heterocycles: {final_product_heterocycles}")
    print(f"- New heterocycles in final product: {heterocycle_difference}")

    condition1 = heterocycle_count >= 3 or (heterocycle_count >= 2 and has_nitrile)
    condition2 = final_product_n_count >= 4
    condition3 = (
        len(synthesized_heterocycles) > 0
        or heterocycle_forming_reactions
        or modified_heterocycles
        or heterocycle_transformation
    )

    print(f"Condition 1 (heterocycle count): {condition1}")
    print(f"Condition 2 (nitrogen count): {condition2}")
    print(f"Condition 3 (synthesis/modification): {condition3}")

    if condition1 and condition2 and condition3:
        found_n_rich_synthesis = True
        print(
            f"Found nitrogen-rich heterocycle synthesis with {heterocycle_count} heterocycles, "
            f"{final_product_n_count} nitrogen atoms"
        )

    # If we have a nitrogen-rich molecule but couldn't detect synthesis, check if the route is complex enough
    # to warrant considering it a synthesis (at least 3 reaction steps)
    if condition1 and condition2 and not condition3:
        # Count reaction steps
        reaction_count = 0

        def count_reactions(node):
            nonlocal reaction_count
            if node["type"] == "reaction":
                reaction_count += 1
            for child in node.get("children", []):
                count_reactions(child)

        count_reactions(route)
        print(f"Route contains {reaction_count} reaction steps")

        # If we have at least 3 reaction steps, consider it a synthesis
        if reaction_count >= 3:
            found_n_rich_synthesis = True
            print(f"Found nitrogen-rich heterocycle synthesis based on route complexity")

    return found_n_rich_synthesis
