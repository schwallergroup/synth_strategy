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
    Detects if the synthesis follows a linear strategy without convergent steps.

    A linear synthesis has each reaction step with only one significant reactant that
    continues the synthetic pathway. Convergent synthesis combines multiple complex
    building blocks in at least one step.

    Args:
        route: A synthesis route JSON object following the SynthesisRoute schema

    Returns:
        bool: True if the synthesis is linear (no convergent steps), False otherwise
    """
    has_convergent_step = False

    # Common reagent-based reactions that shouldn't be considered convergent
    reagent_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Reduction of aldehydes and ketones to alcohols",
        "Boc amine protection",
        "Boc amine deprotection",
        "Esterification of Carboxylic Acids",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Alcohol protection with silyl ethers",
        "Alcohol deprotection from silyl ethers",
        "Methylation",
        "Hydrogenation (double to single)",
        "Hydrogenation (triple to double)",
        "Wittig reaction with triphenylphosphorane",
        "Wittig with Phosphonium",
        "Wittig",
        "Michael addition",
        "Horner-Wadsworth-Emmons",
        "aza-Michael addition primary",
        "aza-Michael addition secondary",
        "aza-Michael addition aromatic",
        "thia-Michael addition",
        "oxa-Michael addition",
    ]

    # Common reagents that shouldn't be considered significant reactants
    reagent_functional_groups = [
        "Triflate",
        "Mesylate",
        "Tosylate",
        "Magnesium halide",
        "Zinc halide",
        "Alkyl lithium",
        "Silyl protective group",
        "Phosphate ester",
        "Boronic acid",
        "Boronic ester",
    ]

    def is_reagent(smiles):
        """Check if a molecule is likely a reagent rather than a building block"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for common reagent functional groups
        if any(checker.check_fg(fg, smiles) for fg in reagent_functional_groups):
            return True

        # Check for small molecules that are common reagents
        atom_count = mol.GetNumAtoms()
        if atom_count <= 6:  # Small molecules are often reagents
            return True

        # Check for simple inorganic compounds
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count == 0:  # Inorganic compounds
            return True

        # Check for phosphorus-containing reagents (often Wittig reagents)
        has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
        if has_phosphorus:
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Check if this is a reagent-based reaction
            is_reagent_reaction = any(
                checker.check_reaction(rxn, rsmi) for rxn in reagent_reactions
            )
            if is_reagent_reaction:
                print(f"Reagent-based reaction detected at depth {depth}: {rsmi}")
                # For reagent-based reactions, we don't consider them convergent
                return

            # Identify significant reactants based on chemical complexity
            significant_reactants = []
            for r in reactants:
                if is_reagent(r):
                    print(f"Reagent detected at depth {depth}: {r}")
                    continue

                mol = Chem.MolFromSmiles(r)
                if mol:
                    # Count carbon atoms (organic molecules typically have carbon backbone)
                    carbon_count = sum(
                        1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6
                    )

                    # Check for presence of functional groups (indicates chemical complexity)
                    has_functional_groups = any(
                        checker.check_fg(fg, r)
                        for fg in [
                            "Carboxylic acid",
                            "Ester",
                            "Amide",
                            "Amine",
                            "Alcohol",
                            "Aldehyde",
                            "Ketone",
                            "Alkene",
                            "Alkyne",
                            "Aromatic halide",
                            "Nitrile",
                            "Nitro group",
                        ]
                    )

                    # Check for presence of ring structures (indicates complexity)
                    has_rings = any(
                        checker.check_ring(ring, r)
                        for ring in [
                            "benzene",
                            "pyridine",
                            "cyclohexane",
                            "cyclopentane",
                            "furan",
                            "thiophene",
                            "pyrrole",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                            "pyrimidine",
                            "piperidine",
                            "morpholine",
                        ]
                    )

                    # Consider a reactant significant if it has sufficient complexity
                    # Either has substantial carbon count OR has both carbon backbone and functional groups/rings
                    if carbon_count >= 7 or (
                        carbon_count >= 4 and (has_functional_groups and has_rings)
                    ):
                        significant_reactants.append(r)
                        print(
                            f"Significant reactant at depth {depth}: {r} (C atoms: {carbon_count}, FG: {has_functional_groups}, Rings: {has_rings})"
                        )

            # If there are multiple significant reactants, this is a convergent step
            # Check for specific reaction types that might appear convergent but are typically linear
            if len(significant_reactants) > 1:
                # Check if this is a Wittig-type or HWE reaction which are often considered linear
                is_wittig_type = any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in [
                        "Wittig reaction with triphenylphosphorane",
                        "Wittig with Phosphonium",
                        "Wittig",
                        "Horner-Wadsworth-Emmons",
                    ]
                )

                if not is_wittig_type:
                    has_convergent_step = True
                    print(
                        f"Found convergent step at depth {depth} with {len(significant_reactants)} significant reactants"
                    )
                else:
                    print(
                        f"Wittig-type reaction at depth {depth} - not considered convergent"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if the synthesis is linear (no convergent steps)
    return not has_convergent_step
