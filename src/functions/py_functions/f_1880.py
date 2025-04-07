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
    Detects if the synthesis route involves nucleophilic aromatic substitution.
    """
    has_snar = False

    def dfs_traverse(node):
        nonlocal has_snar

        if node["type"] == "reaction" and not has_snar:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # First check if this is a known SNAr reaction type
                if (
                    checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                ):
                    print(
                        f"Found nucleophilic aromatic substitution via reaction check: {rsmi}"
                    )
                    has_snar = True
                    return

                # If not a known reaction type, check for characteristic patterns
                reactant_mols = []
                for r in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            reactant_mols.append(mol)
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                except:
                    product_mol = None

                if not product_mol or not reactant_mols:
                    return

                # Check for aromatic halides in reactants
                has_aromatic_halide = False
                for r_smiles in reactants_smiles:
                    if checker.check_fg("Aromatic halide", r_smiles):
                        has_aromatic_halide = True
                        break

                # Check for activating groups in reactants
                has_activating_group = False
                for r_smiles in reactants_smiles:
                    if (
                        checker.check_fg("Nitro group", r_smiles)
                        or checker.check_fg("Nitrile", r_smiles)
                        or checker.check_fg("Ester", r_smiles)
                        or checker.check_fg("Ketone", r_smiles)
                        or checker.check_ring("pyridine", r_smiles)
                        or checker.check_ring("pyrimidine", r_smiles)
                        or checker.check_ring("pyrazine", r_smiles)
                    ):
                        has_activating_group = True
                        break

                # Check for nucleophiles in reactants
                has_nucleophile = False
                for r_smiles in reactants_smiles:
                    if (
                        checker.check_fg("Primary amine", r_smiles)
                        or checker.check_fg("Secondary amine", r_smiles)
                        or checker.check_fg("Phenol", r_smiles)
                        or checker.check_fg("Primary alcohol", r_smiles)
                        or checker.check_fg("Secondary alcohol", r_smiles)
                        or checker.check_fg("Tertiary alcohol", r_smiles)
                        or checker.check_fg("Aromatic thiol", r_smiles)
                        or checker.check_fg("Aliphatic thiol", r_smiles)
                    ):
                        has_nucleophile = True
                        break

                # Check if the product has a new C-N, C-O, or C-S bond
                # This is a simplification - ideally we would use atom mapping to track the exact substitution
                has_new_bond = False
                for r_smiles in reactants_smiles:
                    if checker.check_fg(
                        "Aromatic halide", r_smiles
                    ) and not checker.check_fg("Aromatic halide", product_smiles):
                        has_new_bond = True
                        break

                if (
                    has_aromatic_halide
                    and has_activating_group
                    and has_nucleophile
                    and has_new_bond
                ):
                    print(
                        f"Found nucleophilic aromatic substitution via pattern matching: {rsmi}"
                    )
                    has_snar = True

            except Exception as e:
                print(f"Error in SNAr detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar
