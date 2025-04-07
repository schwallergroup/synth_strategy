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


def main(route, threshold=3):
    """
    Detects if the synthesis involves multiple C-N bond formations.

    Args:
        route: The synthesis route to analyze
        threshold: Minimum number of C-N bond formations required (default: 3)

    Returns:
        bool: True if the synthesis involves at least 'threshold' C-N bond formations
    """
    cn_bond_formations = 0

    def dfs_traverse(node):
        nonlocal cn_bond_formations

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for common C-N bond formation reactions
                cn_bond_formed = False

                # Check for amide formation reactions
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with ammonia to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    )
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Ester with secondary amine to amide", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Amide formation")

                # Check for N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)
                elif (
                    checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        rsmi,
                    )
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "Goldberg coupling aryl amine-aryl chloride", rsmi
                    )
                    or checker.check_reaction(
                        "Goldberg coupling aryl amide-aryl chloride", rsmi
                    )
                    or checker.check_reaction("Goldberg coupling", rsmi)
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution amine", rsmi
                    )
                    or checker.check_reaction("N-arylation_heterocycles", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: N-arylation")

                # Check for reductive amination
                elif (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Reductive amination")

                # Check for alkylation of amines
                elif (
                    checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("N-methylation", rsmi)
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                    or checker.check_reaction(
                        "Eschweiler-Clarke Primary Amine Methylation", rsmi
                    )
                    or checker.check_reaction(
                        "Eschweiler-Clarke Secondary Amine Methylation", rsmi
                    )
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Amine alkylation")

                # Check for urea/thiourea formation
                elif (
                    checker.check_reaction(
                        "Urea synthesis via isocyanate and primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and diazo", rsmi
                    )
                    or checker.check_reaction(
                        "Urea synthesis via isocyanate and sulfonamide", rsmi
                    )
                    or checker.check_reaction("urea", rsmi)
                    or checker.check_reaction("thiourea", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Urea/thiourea formation")

                # Check for heterocycle formation involving C-N bonds
                elif (
                    checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                    or checker.check_reaction("Paal-Knorr pyrrole synthesis", rsmi)
                    or checker.check_reaction("Paal-Knorr pyrrole", rsmi)
                    or checker.check_reaction("Pictet-Spengler", rsmi)
                    or checker.check_reaction(
                        "benzimidazole_derivatives_carboxylic-acid/ester", rsmi
                    )
                    or checker.check_reaction(
                        "benzimidazole_derivatives_aldehyde", rsmi
                    )
                    or checker.check_reaction(
                        "Benzimidazole formation from aldehyde", rsmi
                    )
                    or checker.check_reaction(
                        "Benzimidazole formation from acyl halide", rsmi
                    )
                    or checker.check_reaction(
                        "Benzimidazole formation from ester/carboxylic acid", rsmi
                    )
                    or checker.check_reaction("Benzimidazole aldehyde", rsmi)
                    or checker.check_reaction("imidazole", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Heterocycle formation")

                # Check for imine formation
                elif (
                    checker.check_reaction(
                        "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of primary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to ketones/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction(
                        "Addition of secondary amines to aldehydes/thiocarbonyls", rsmi
                    )
                    or checker.check_reaction("Ketone/aldehyde to hydrazone", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Imine formation")

                # Check for sulfonamide formation
                elif (
                    checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )
                    or checker.check_reaction("sulfon_amide", rsmi)
                ):
                    cn_bond_formed = True
                    print(f"Found C-N bond formation: Sulfonamide formation")

                # If none of the specific reactions were detected, try a more general approach
                if not cn_bond_formed:
                    # Create molecules for reactants and product
                    reactants = reactants_part.split(".")
                    reactant_mols = [
                        Chem.MolFromSmiles(r)
                        for r in reactants
                        if r and Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and all(reactant_mols):
                        # Check for amine functional groups in reactants
                        amine_in_reactants = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants
                        )

                        # Check for carbon-containing functional groups that might form C-N bonds
                        carbon_fg_in_reactants = any(
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            or checker.check_fg("Aromatic halide", r)
                            or checker.check_fg("Carboxylic acid", r)
                            or checker.check_fg("Ester", r)
                            or checker.check_fg("Acyl halide", r)
                            or checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Ketone", r)
                            for r in reactants
                        )

                        # Check for amide, amine, or other C-N containing groups in product
                        cn_in_product = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                            or checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                            or checker.check_fg("Urea", product)
                            or checker.check_fg("Thiourea", product)
                            or checker.check_fg("Sulfonamide", product)
                        )

                        if (
                            amine_in_reactants
                            and carbon_fg_in_reactants
                            and cn_in_product
                        ):
                            cn_bond_formed = True
                            print(f"Found C-N bond formation: General pattern")

                if cn_bond_formed:
                    cn_bond_formations += 1

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total C-N bond formations found: {cn_bond_formations}")
    return cn_bond_formations >= threshold
