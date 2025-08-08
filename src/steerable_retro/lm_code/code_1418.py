#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    Detects late-stage amine alkylation via reductive amination or other alkylation methods in the final steps.
    """
    late_stage_alkylations = 0

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_alkylations

        if node["type"] == "reaction" and depth <= 2:  # Focus on late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for patterns consistent with amine alkylation
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants if r)
                has_formaldehyde = any(checker.check_fg("Formaldehyde", r) for r in reactants if r)
                # Additional check for common formaldehyde SMILES patterns
                if not has_formaldehyde:
                    has_formaldehyde = any(
                        r.strip() in ["O=CH2", "C=O", "CH2O"] for r in reactants if r
                    )

                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants if r)
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants if r
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants if r
                )
                has_alkyl_halide = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    for r in reactants
                    if r
                )

                print(f"  Has aldehyde: {has_aldehyde}")
                print(f"  Has formaldehyde: {has_formaldehyde}")
                print(f"  Has ketone: {has_ketone}")
                print(f"  Has primary amine: {has_primary_amine}")
                print(f"  Has secondary amine: {has_secondary_amine}")
                print(f"  Has alkyl halide: {has_alkyl_halide}")

                # Check if product has more substituted amine than reactants
                product_has_secondary_amine = checker.check_fg("Secondary amine", product)
                product_has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                print(f"  Product has secondary amine: {product_has_secondary_amine}")
                print(f"  Product has tertiary amine: {product_has_tertiary_amine}")

                # Check for various amine alkylation reactions
                is_reductive_amination_aldehyde = checker.check_reaction(
                    "Reductive amination with aldehyde", rsmi
                )
                is_reductive_amination_ketone = checker.check_reaction(
                    "Reductive amination with ketone", rsmi
                )
                is_reductive_amination_alcohol = checker.check_reaction(
                    "Reductive amination with alcohol", rsmi
                )
                is_n_alkylation = checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rsmi
                ) or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                is_methylation = (
                    checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                    or checker.check_reaction("DMS Amine methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                    or checker.check_reaction("N-methylation", rsmi)
                )

                # Check for reductive amination with formaldehyde specifically
                is_reductive_amination_formaldehyde = (
                    has_formaldehyde
                    and (has_primary_amine or has_secondary_amine)
                    and (product_has_secondary_amine or product_has_tertiary_amine)
                )

                print(f"  Is reductive amination with aldehyde: {is_reductive_amination_aldehyde}")
                print(f"  Is reductive amination with ketone: {is_reductive_amination_ketone}")
                print(f"  Is reductive amination with alcohol: {is_reductive_amination_alcohol}")
                print(
                    f"  Is reductive amination with formaldehyde: {is_reductive_amination_formaldehyde}"
                )
                print(f"  Is N-alkylation: {is_n_alkylation}")
                print(f"  Is methylation: {is_methylation}")

                # Detect amine alkylation
                is_amine_alkylation = False

                # Case 1: Reductive amination
                if (
                    (has_aldehyde or has_ketone or has_formaldehyde)
                    and (has_primary_amine or has_secondary_amine)
                    and (
                        is_reductive_amination_aldehyde
                        or is_reductive_amination_ketone
                        or is_reductive_amination_alcohol
                        or is_reductive_amination_formaldehyde
                    )
                ):
                    is_amine_alkylation = True

                # Case 2: Direct N-alkylation
                elif (
                    (has_primary_amine or has_secondary_amine)
                    and (has_alkyl_halide or has_formaldehyde)
                    and (is_n_alkylation or is_methylation)
                ):
                    is_amine_alkylation = True

                # Case 3: Any methylation reaction that produces a more substituted amine
                elif is_methylation and (
                    (has_primary_amine and product_has_secondary_amine)
                    or (has_secondary_amine and product_has_tertiary_amine)
                ):
                    is_amine_alkylation = True

                # Case 4: Pattern-based detection for formaldehyde methylation
                elif has_formaldehyde and (
                    (has_primary_amine and product_has_secondary_amine)
                    or (has_secondary_amine and product_has_tertiary_amine)
                ):
                    # Verify that a methyl group is added (formaldehyde -> methyl)
                    is_amine_alkylation = True
                    print(f"  Detected formaldehyde methylation pattern")

                if is_amine_alkylation:
                    late_stage_alkylations += 1
                    print(f"Found late-stage amine alkylation at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Total late-stage amine alkylations found: {late_stage_alkylations}")

    return late_stage_alkylations >= 1
