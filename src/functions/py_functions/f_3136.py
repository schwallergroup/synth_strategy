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
    Detects if the synthesis proceeds without any ring formation or opening steps
    """
    has_ring_change = False

    # List of common ring types to check
    ring_types = [
        "furan",
        "pyran",
        "dioxane",
        "tetrahydrofuran",
        "tetrahydropyran",
        "oxirane",
        "oxetane",
        "oxolane",
        "oxane",
        "dioxolane",
        "dioxolene",
        "trioxane",
        "dioxepane",
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
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
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
        "thiophene",
        "thiopyran",
        "thiirane",
        "thietane",
        "thiolane",
        "thiane",
        "dithiane",
        "dithiolane",
        "benzothiophene",
        "oxathiolane",
        "dioxathiolane",
        "thiazolidine",
        "oxazolidine",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "cycloheptane",
        "cyclooctane",
        "benzene",
        "naphthalene",
        "anthracene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "pteridin",
        "phenothiazine",
        "phenoxazine",
        "dibenzofuran",
        "dibenzothiophene",
        "xanthene",
        "thioxanthene",
        "pyrroline",
        "pyrrolidone",
        "imidazolidine",
        "porphyrin",
        "indazole",
        "benzotriazole",
    ]

    # List of known ring-forming and ring-opening reactions
    ring_changing_reactions = [
        # Ring-forming reactions
        "Paal-Knorr pyrrole synthesis",
        "Formation of NOS Heterocycles",
        "Diels-Alder",
        "Pauson-Khand reaction",
        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
        "Huisgen 1,3 dipolar cycloaddition",
        "Huisgen alkene-azide 1,3 dipolar cycloaddition",
        "Pyrazole formation",
        "Intramolecular amination of azidobiphenyls (heterocycle formation)",
        "Intramolecular amination (heterocycle formation)",
        "Pictet-Spengler",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "spiro-chromanone",
        "pyrazole",
        "phthalazinone",
        "Paal-Knorr pyrrole",
        "triaryl-imidazole",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "piperidine_indole",
        "imidazole",
        # Ring-opening reactions
        "Acetal hydrolysis to diol",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Ring opening of epoxide with amine",
        "Retro-Diels-Alder from oxazole",
    ]

    # Reactions that may change ring count but aren't considered ring-forming/breaking
    exempt_reactions = [
        "Hydrogenolysis of amides/imides/carbamates",
        "Hydrogenolysis of tertiary amines",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "Hydrogenolysis",
    ]

    def dfs_traverse(node):
        nonlocal has_ring_change

        # Skip further processing if we already found a ring change
        if has_ring_change:
            return

        # Process the current node
        if node["type"] == "reaction" and "metadata" in node:
            # Check if reaction is explicitly marked as ring-breaking
            if node["metadata"].get("RingBreaker", False):
                print(f"Found reaction explicitly marked as RingBreaker in metadata")
                has_ring_change = True
                return

            # Check reaction SMILES if available
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for exempt reactions first
                for reaction_type in exempt_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Skipping exempt reaction: {reaction_type} in reaction: {rsmi}"
                        )
                        # Continue with traversal but don't mark as ring change
                        break
                else:  # This else belongs to the for loop - executes if no break occurred
                    # Check for known ring-changing reaction types
                    for reaction_type in ring_changing_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Found ring-changing reaction: {reaction_type} in reaction: {rsmi}"
                            )
                            has_ring_change = True
                            return

                    # If not a known reaction type, analyze the reaction components
                    reactants_part = rsmi.split(">")[0]
                    product = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check for hydrogenolysis pattern (common in benzyl deprotection)
                    if "[Pd+2]" in rsmi or "Pd" in rsmi or "H2" in rsmi:
                        # Check if this is a benzyl deprotection
                        benzyl_pattern = False
                        for reactant in reactants:
                            if checker.check_ring("benzene", reactant):
                                # Check if benzene is connected to a functional group that's being removed
                                if (
                                    checker.check_fg("Primary amine", product)
                                    or checker.check_fg("Secondary amine", product)
                                    or checker.check_fg("Primary alcohol", product)
                                    or checker.check_fg("Secondary alcohol", product)
                                    or checker.check_fg("Carboxylic acid", product)
                                ):
                                    benzyl_pattern = True
                                    break

                        if benzyl_pattern:
                            print(f"Skipping benzyl deprotection reaction: {rsmi}")
                            # Continue with traversal but don't mark as ring change
                            for child in node.get("children", []):
                                dfs_traverse(child)
                            return

                    # Check specific ring types in reactants and products
                    reactant_rings = set()
                    for reactant in reactants:
                        for ring_type in ring_types:
                            if checker.check_ring(ring_type, reactant):
                                reactant_rings.add(ring_type)

                    product_rings = set()
                    for ring_type in ring_types:
                        if checker.check_ring(ring_type, product):
                            product_rings.add(ring_type)

                    # If different ring types are present, it's a ring transformation
                    if reactant_rings != product_rings:
                        # Check if the only difference is benzene rings (could be benzyl deprotection)
                        diff_reactant = reactant_rings - product_rings
                        diff_product = product_rings - reactant_rings

                        if diff_reactant == {"benzene"} and not diff_product:
                            # This might be a benzyl deprotection - check if it's a protection/deprotection reaction
                            if any(
                                checker.check_reaction(rxn, rsmi)
                                for rxn in [
                                    "Hydroxyl benzyl deprotection",
                                    "Carboxyl benzyl deprotection",
                                    "Boc amine deprotection",
                                ]
                            ):
                                print(f"Skipping benzyl deprotection: {rsmi}")
                                # Continue with traversal but don't mark as ring change
                                for child in node.get("children", []):
                                    dfs_traverse(child)
                                return

                        print(
                            f"Found ring structure change: {reactant_rings} -> {product_rings} in reaction: {rsmi}"
                        )
                        has_ring_change = True
                        return

                    # As a fallback, check total ring count
                    reactant_rings_count = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            reactant_rings_count += mol.GetRingInfo().NumRings()

                    product_mol = Chem.MolFromSmiles(product)
                    product_rings_count = 0
                    if product_mol:
                        product_rings_count = product_mol.GetRingInfo().NumRings()

                    # If ring count changes significantly (more than just benzyl removal)
                    if (
                        abs(reactant_rings_count - product_rings_count) > 0
                        and abs(reactant_rings_count - product_rings_count) != 1
                    ):
                        print(
                            f"Found significant ring count change: {reactant_rings_count} -> {product_rings_count} in reaction: {rsmi}"
                        )
                        has_ring_change = True
                        return

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if no ring changes were detected
    return not has_ring_change
