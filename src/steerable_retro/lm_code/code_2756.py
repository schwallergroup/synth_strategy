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
    Detects if the synthetic route involves late-stage N-methylation via reductive amination.
    """
    n_methylation = False

    def dfs_traverse(node, depth=0):
        nonlocal n_methylation

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Late stage is defined as depth ≤ 2
            if depth <= 2:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Direct check for N-methylation reactions
                is_n_methylation = (
                    checker.check_reaction("N-methylation", rsmi)
                    or checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("DMS Amine methylation", rsmi)
                )

                # Check for reductive amination that results in N-methylation
                if not is_n_methylation:
                    # Check for reductive amination reactions
                    is_reductive_amination = (
                        checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction("reductive amination", rsmi)
                        or checker.check_reaction("{reductive amination}", rsmi)
                    )

                    # Check for formaldehyde in reactants
                    has_formaldehyde = (
                        any(checker.check_fg("Formaldehyde", r) for r in reactants)
                        or "O=[CH2]" in rsmi.split(">")[0]
                        or "C=O" in rsmi.split(">")[0]
                    )

                    # Check for amine in reactants
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )

                    # Check for N-methyl in product
                    has_n_methyl_product = checker.check_fg(
                        "Tertiary amine", product
                    ) or checker.check_fg("Secondary amine", product)

                    # Check for N-CH3 formation by analyzing atom mapping
                    # Look for patterns like [CH3:x][N:y] in product where [CH:x] and [N:y] are separate in reactants
                    reactants_str = rsmi.split(">")[0]
                    product_str = rsmi.split(">")[-1]

                    # Check if there's a formaldehyde carbon ([CH2]) that becomes attached to a nitrogen
                    import re

                    formaldehyde_carbons = re.findall(r"\[CH2:(\d+)\]", reactants_str)
                    nitrogen_atoms = re.findall(r"\[N[H]?:(\d+)\]", reactants_str)

                    # Check if any formaldehyde carbon is attached to any nitrogen in the product
                    n_methylation_pattern = False
                    for c_idx in formaldehyde_carbons:
                        for n_idx in nitrogen_atoms:
                            if (
                                f"[CH3:{c_idx}][N:{n_idx}]" in product_str
                                or f"[CH3:{c_idx}][N" in product_str
                            ):
                                n_methylation_pattern = True
                                print(
                                    f"Found N-methylation pattern: Carbon {c_idx} attached to Nitrogen {n_idx}"
                                )
                                break

                    # If we have a reductive amination with formaldehyde and an amine, and the product has N-methyl
                    # or we detected a direct N-methylation pattern, mark as N-methylation
                    if (
                        is_reductive_amination
                        and has_formaldehyde
                        and has_amine
                        and has_n_methyl_product
                    ) or n_methylation_pattern:
                        is_n_methylation = True

                if is_n_methylation:
                    print(
                        f"Detected late-stage N-methylation at depth {depth}. Reactants: {reactants}, Product: {product}"
                    )
                    n_methylation = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_methylation
