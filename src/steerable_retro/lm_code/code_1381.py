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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a strategy involving early nitration followed by nucleophilic aromatic substitution.
    """
    # Initialize tracking variables
    nitration_depth = -1
    chlorination_depth = -1
    substitution_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitration_depth, chlorination_depth, substitution_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitration
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                ):
                    nitration_depth = depth
                    print(f"Detected nitration at depth {depth}, product: {product}")

                # Check for chlorination
                if checker.check_reaction("Aromatic chlorination", rsmi):
                    chlorination_depth = depth
                    print(f"Detected chlorination at depth {depth}, product: {product}")

                # Check for nucleophilic aromatic substitution
                if (
                    checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                ):
                    substitution_depth = depth
                    print(
                        f"Detected nucleophilic aromatic substitution at depth {depth}, product: {product}"
                    )

                # If we don't have specific reaction checkers, fall back to functional group analysis
                if (
                    nitration_depth == -1
                    and checker.check_fg("Nitro group", product)
                    and not any(checker.check_fg("Nitro group", r) for r in reactants)
                ):
                    nitration_depth = depth
                    print(f"Detected nitration (by FG) at depth {depth}, product: {product}")

                if (
                    chlorination_depth == -1
                    and checker.check_fg("Aromatic halide", product)
                    and not any(checker.check_fg("Aromatic halide", r) for r in reactants)
                ):
                    # Verify it's specifically chlorination
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c-Cl")):
                        chlorination_depth = depth
                        print(f"Detected chlorination (by FG) at depth {depth}, product: {product}")

                if substitution_depth == -1:
                    # Check for C-Cl to C-N transformation
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product)

                    if all(reactant_mols) and product_mol:
                        has_cl_reactant = any(
                            mol.HasSubstructMatch(Chem.MolFromSmarts("c-Cl"))
                            for mol in reactant_mols
                        )
                        has_nitro_reactant = any(
                            checker.check_fg("Nitro group", Chem.MolToSmiles(mol))
                            for mol in reactant_mols
                        )
                        has_cl_product = product_mol.HasSubstructMatch(Chem.MolFromSmarts("c-Cl"))
                        has_nh_product = product_mol.HasSubstructMatch(Chem.MolFromSmarts("c-[N]"))

                        if (
                            has_cl_reactant
                            and has_nitro_reactant
                            and not has_cl_product
                            and has_nh_product
                        ):
                            substitution_depth = depth
                            print(
                                f"Detected nucleophilic aromatic substitution (by FG) at depth {depth}, product: {product}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Nitration depth: {nitration_depth}, Chlorination depth: {chlorination_depth}, Substitution depth: {substitution_depth}"
    )

    # Check if all three features are detected in the correct sequence
    if nitration_depth != -1 and chlorination_depth != -1 and substitution_depth != -1:
        # Verify correct sequence: nitration (early) -> chlorination -> nucleophilic substitution (late)
        if nitration_depth > chlorination_depth > substitution_depth:
            print("Detected early nitration with nucleophilic aromatic substitution strategy")
            return True

    return False
