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
    Detects if the synthesis route includes oxidation of an aromatic methyl group to a carboxylic acid.
    """
    methyl_oxidation_found = False

    def dfs_traverse(node):
        nonlocal methyl_oxidation_found

        if node["type"] == "reaction" and not methyl_oxidation_found:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an oxidation reaction
                is_oxidation = (
                    checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of aldehyde to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                    or checker.check_reaction("Aryl halide to carboxylic acid", rsmi)
                )

                if is_oxidation:
                    # Check for aromatic methyl in reactants and carboxylic acid in product
                    product_mol = Chem.MolFromSmiles(product)

                    if product_mol and checker.check_fg("Carboxylic acid", product):
                        print(f"Found product with carboxylic acid: {product}")

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)

                            # Check for aromatic methyl group in reactant
                            aromatic_methyl_pattern = Chem.MolFromSmarts("c-[CH3]")
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                aromatic_methyl_pattern
                            ):
                                print(f"Found reactant with aromatic methyl: {reactant}")

                                # Try to match atom mapping between reactant methyl and product acid
                                # Extract atom-mapped SMILES for analysis
                                reactant_matches = reactant_mol.GetSubstructMatches(
                                    aromatic_methyl_pattern
                                )

                                # Check if the number of methyl groups decreased and carboxylic acids increased
                                methyl_count_reactant = len(
                                    reactant_mol.GetSubstructMatches(aromatic_methyl_pattern)
                                )

                                # Count methyl groups in product
                                product_methyl_count = 0
                                if product_mol.HasSubstructMatch(aromatic_methyl_pattern):
                                    product_methyl_count = len(
                                        product_mol.GetSubstructMatches(aromatic_methyl_pattern)
                                    )

                                # If there's one less methyl in product, likely our transformation
                                if methyl_count_reactant > product_methyl_count:
                                    print("Methyl count decreased from reactant to product")

                                    # Check if the carboxylic acid is attached to an aromatic ring
                                    aromatic_acid_pattern = Chem.MolFromSmarts("c-C(=O)O")
                                    if product_mol.HasSubstructMatch(aromatic_acid_pattern):
                                        print("Found aromatic carboxylic acid in product")
                                        methyl_oxidation_found = True
                                        break

                # Also check for direct methyl to acid oxidation
                # This is a more specific check for the exact transformation we're looking for
                for reactant in reactants:
                    if checker.check_fg("Aromatic halide", reactant) or not checker.check_fg(
                        "Carboxylic acid", reactant
                    ):
                        if checker.check_fg("Carboxylic acid", product):
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                aromatic_methyl_pattern = Chem.MolFromSmarts("c-[CH3]")
                                if reactant_mol.HasSubstructMatch(aromatic_methyl_pattern):
                                    print(f"Direct methyl to acid oxidation detected: {rsmi}")
                                    methyl_oxidation_found = True
                                    break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return methyl_oxidation_found
