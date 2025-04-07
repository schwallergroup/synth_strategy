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
    This function detects if the synthesis route involves heterocycle formation followed by
    sequential functionalization (chlorination, amination).
    """
    # Lists to store detected reactions with their depths
    heterocycle_formations = []
    chlorinations = []
    aminations = []

    # List of heterocycles to check
    heterocycle_types = [
        "pyridine",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "furan",
        "thiophene",
        "morpholine",
        "piperidine",
        "piperazine",
        "pyrrolidine",
    ]

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for heterocycle formation
                product_has_heterocycle = False
                reactants_have_heterocycle = False

                # Check if product contains a heterocycle
                for heterocycle in heterocycle_types:
                    if checker.check_ring(heterocycle, product_smiles):
                        product_has_heterocycle = True
                        print(f"Found {heterocycle} in product at depth {depth}")
                        break

                # Check if any reactant contains a heterocycle
                for reactant in reactants_smiles:
                    for heterocycle in heterocycle_types:
                        if checker.check_ring(heterocycle, reactant):
                            reactants_have_heterocycle = True
                            print(f"Found {heterocycle} in reactant at depth {depth}")
                            break
                    if reactants_have_heterocycle:
                        break

                # Detect heterocycle formation
                if product_has_heterocycle and not reactants_have_heterocycle:
                    heterocycle_formations.append(depth)
                    print(
                        f"Detected heterocycle formation at depth {depth}, product: {product_smiles}"
                    )

                # Also check for specific heterocycle formation reactions
                heterocycle_formation_reactions = [
                    "Formation of NOS Heterocycles",
                    "Paal-Knorr pyrrole synthesis",
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "benzothiazole",
                    "benzoxazole_arom-aldehyde",
                    "benzoxazole_carboxylic-acid",
                    "thiazole",
                    "Benzimidazole formation from aldehyde",
                    "Benzimidazole formation from acyl halide",
                    "Benzimidazole formation from ester/carboxylic acid",
                    "Benzoxazole formation from aldehyde",
                    "Benzoxazole formation from acyl halide",
                    "Benzoxazole formation from ester/carboxylic acid",
                    "Benzoxazole formation (intramolecular)",
                    "Benzothiazole formation from aldehyde",
                    "Benzothiazole formation from acyl halide",
                    "Benzothiazole formation from ester/carboxylic acid",
                ]

                for reaction_name in heterocycle_formation_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        heterocycle_formations.append(depth)
                        print(
                            f"Detected heterocycle formation reaction {reaction_name} at depth {depth}"
                        )
                        break

                # Check for chlorination reaction
                chlorination_reactions = [
                    "Aromatic chlorination",
                    "Chlorination",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                    "Primary amine to chloride",
                ]

                for reaction_name in chlorination_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        # Verify product has chlorine-containing functional group
                        if (
                            checker.check_fg("Aromatic halide", product_smiles)
                            or checker.check_fg("Primary halide", product_smiles)
                            or checker.check_fg("Secondary halide", product_smiles)
                            or checker.check_fg("Tertiary halide", product_smiles)
                        ):
                            chlorinations.append(depth)
                            print(
                                f"Detected chlorination reaction {reaction_name} at depth {depth}"
                            )
                            break

                # Check for chlorinating agents in reactants
                chlorinating_agents = ["POCl3", "SOCl2", "PCl5", "PCl3", "Cl2", "NCS", "CCl4"]
                for reactant in reactants_smiles:
                    if (
                        any(agent in reactant for agent in chlorinating_agents)
                        or "O=P(Cl)(Cl)Cl" in reactant
                    ):
                        # Check if product contains chlorine
                        if "Cl" in product_smiles and (
                            checker.check_fg("Aromatic halide", product_smiles)
                            or checker.check_fg("Primary halide", product_smiles)
                            or checker.check_fg("Secondary halide", product_smiles)
                            or checker.check_fg("Tertiary halide", product_smiles)
                        ):
                            chlorinations.append(depth)
                            print(f"Detected chlorination with agent at depth {depth}")
                            break

                # Check for direct chlorine addition by comparing reactants and products
                if not chlorinations or depth not in chlorinations:
                    # Check if product has chlorine and reactant has OH or =O that could be replaced
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and "Cl" in product_smiles:
                        # Check for pyridine chlorination specifically
                        if checker.check_ring("pyridine", product_smiles) and any(
                            "pyridine" in r for r in reactants_smiles
                        ):
                            for reactant in reactants_smiles:
                                if "=O" in reactant and checker.check_ring("pyridine", reactant):
                                    chlorinations.append(depth)
                                    print(f"Detected pyridine chlorination at depth {depth}")
                                    break

                        # General check for chlorination
                        if (
                            checker.check_fg("Aromatic halide", product_smiles)
                            or checker.check_fg("Primary halide", product_smiles)
                            or checker.check_fg("Secondary halide", product_smiles)
                            or checker.check_fg("Tertiary halide", product_smiles)
                        ):
                            chlorinations.append(depth)
                            print(
                                f"Detected chlorine addition at depth {depth} by comparing reactants and products"
                            )

                # Check for amination reaction
                amination_reactions = [
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Alkylation of amines",
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Reduction of nitro groups to amines",
                    "Reduction of nitrile to amine",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Azide to amine reduction (Staudinger)",
                    "Buchwald-Hartwig",
                    "Goldberg coupling aryl amine-aryl chloride",
                    "Ullmann-Goldberg Substitution amine",
                ]

                for reaction_name in amination_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        # Verify product has amine group
                        if (
                            checker.check_fg("Primary amine", product_smiles)
                            or checker.check_fg("Secondary amine", product_smiles)
                            or checker.check_fg("Tertiary amine", product_smiles)
                            or checker.check_fg("Aniline", product_smiles)
                        ):
                            aminations.append(depth)
                            print(f"Detected amination reaction {reaction_name} at depth {depth}")
                            break

                # Check for direct halide to amine transformation
                if not aminations or depth not in aminations:
                    for reactant in reactants_smiles:
                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ) and (
                            checker.check_fg("Primary amine", product_smiles)
                            or checker.check_fg("Secondary amine", product_smiles)
                            or checker.check_fg("Tertiary amine", product_smiles)
                            or checker.check_fg("Aniline", product_smiles)
                        ):
                            aminations.append(depth)
                            print(f"Detected halide to amine transformation at depth {depth}")
                            break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if we found all three reaction types
    if not (heterocycle_formations and chlorinations and aminations):
        print("Missing one or more required reaction types:")
        print(f"Heterocycle formations: {heterocycle_formations}")
        print(f"Chlorinations: {chlorinations}")
        print(f"Aminations: {aminations}")
        return False

    # In retrosynthesis, later stages have lower depth values
    # So we need to find min depth for each reaction type
    min_heterocycle_depth = min(heterocycle_formations)
    min_chlorination_depth = min(chlorinations)
    min_amination_depth = min(aminations)

    # Check if the sequence is in the correct order (in retrosynthesis)
    # Amination (late stage) -> Chlorination -> Heterocycle formation (early stage)
    # So in terms of depth: heterocycle_depth > chlorination_depth > amination_depth
    correct_sequence = min_heterocycle_depth > min_chlorination_depth > min_amination_depth

    print(f"Heterocycle formation depth: {min_heterocycle_depth}")
    print(f"Chlorination depth: {min_chlorination_depth}")
    print(f"Amination depth: {min_amination_depth}")
    print(f"Correct sequence: {correct_sequence}")

    return correct_sequence
