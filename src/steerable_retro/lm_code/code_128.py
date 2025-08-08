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
    Detects a synthetic strategy involving halogenation of an aromatic ring,
    followed by formylation (replacing halogen with aldehyde),
    followed by condensation to form a C=N bond.
    """
    # Track if we found each step in the sequence
    halogenation_depth = None
    formylation_depth = None
    condensation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal halogenation_depth, formylation_depth, condensation_depth

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES at depth {depth}")
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Processing reaction at depth {depth}: {rsmi}")

                # Check for halogenation (aromatic ring gets halogenated)
                if (
                    checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                    or (
                        any(not checker.check_fg("Aromatic halide", r) for r in reactants_smiles)
                        and checker.check_fg("Aromatic halide", product_smiles)
                    )
                ):
                    print(f"Found halogenation reaction at depth {depth}")
                    halogenation_depth = depth

                # Check for formylation (replacing halogen with aldehyde)
                if (
                    any(checker.check_fg("Aromatic halide", r) for r in reactants_smiles)
                    and checker.check_fg("Aldehyde", product_smiles)
                    and not any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)
                ):
                    print(f"Found formylation reaction at depth {depth}")
                    formylation_depth = depth

                # Check for condensation (aldehyde + amine/hydrazine → C=N bond)
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)
                has_hydrazine = any(checker.check_fg("Hydrazine", r) for r in reactants_smiles)
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants_smiles
                )

                if has_aldehyde and (has_hydrazine or has_primary_amine):
                    # Check if product has hydrazone or imine (C=N bond)
                    if (
                        checker.check_fg("Hydrazone", product_smiles)
                        or checker.check_fg("Substituted imine", product_smiles)
                        or checker.check_fg("Unsubstituted imine", product_smiles)
                        or checker.check_reaction(
                            "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                        )
                        or checker.check_reaction("Ketone/aldehyde to hydrazone", rsmi)
                    ):
                        print(f"Found condensation reaction at depth {depth}")
                        condensation_depth = depth

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Halogenation depth: {halogenation_depth}")
    print(f"Formylation depth: {formylation_depth}")
    print(f"Condensation depth: {condensation_depth}")

    # Check if all steps are found and in the correct sequence
    # In retrosynthetic traversal, lower depth means later stage
    if (
        halogenation_depth is not None
        and formylation_depth is not None
        and condensation_depth is not None
        and halogenation_depth > formylation_depth > condensation_depth
    ):
        print("All steps found in correct sequence")
        return True

    return False
