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
    This function detects a synthetic strategy involving propargylation of a ketone
    to form a propargyl alcohol.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    has_propargylation = False

    def dfs_traverse(node):
        nonlocal has_propargylation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for ketone in reactants
                    ketone_in_reactants = False
                    ketone_reactant = None
                    for reactant in reactants:
                        if checker.check_fg("Ketone", reactant):
                            ketone_in_reactants = True
                            ketone_reactant = reactant
                            break

                    # Check for propargyl alcohol in product
                    has_alcohol = checker.check_fg(
                        "Secondary alcohol", product
                    ) or checker.check_fg("Tertiary alcohol", product)
                    has_alkyne = checker.check_fg("Alkyne", product)

                    # Check for propargyl reagent in reactants
                    has_propargyl_reagent = False
                    for reactant in reactants:
                        if checker.check_fg("Alkyne", reactant):
                            # Check if it's a terminal alkyne or propargyl-like structure
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # Check for terminal alkyne or propargyl group
                                terminal_alkyne_pattern = Chem.MolFromSmarts("C#C[H]")
                                propargyl_pattern = Chem.MolFromSmarts("C#CC")
                                if mol.HasSubstructMatch(
                                    terminal_alkyne_pattern
                                ) or mol.HasSubstructMatch(propargyl_pattern):
                                    has_propargyl_reagent = True
                                    break

                    # Check for relevant reaction types
                    is_relevant_reaction = (
                        checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                        or checker.check_reaction(
                            "Reaction of alkyl halides with organometallic coumpounds", rsmi
                        )
                        or checker.check_reaction(
                            "Addition of primary amines to ketones/thiocarbonyls", rsmi
                        )
                        or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                        or checker.check_reaction("Grignard with CO2 to carboxylic acid", rsmi)
                        or checker.check_reaction(
                            "Olefination of ketones with Grignard reagents", rsmi
                        )
                    )

                    # Check if the alkyne is near the alcohol in the product
                    is_propargyl_alcohol = False
                    if has_alcohol and has_alkyne:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Look for propargyl alcohol pattern (alkyne within 3 bonds of OH)
                            propargyl_alcohol_pattern = Chem.MolFromSmarts("C#CC[OH]")
                            propargyl_alcohol_pattern2 = Chem.MolFromSmarts("C#CC([OH])")
                            propargyl_alcohol_pattern3 = Chem.MolFromSmarts("C#CCC[OH]")
                            if (
                                product_mol.HasSubstructMatch(propargyl_alcohol_pattern)
                                or product_mol.HasSubstructMatch(propargyl_alcohol_pattern2)
                                or product_mol.HasSubstructMatch(propargyl_alcohol_pattern3)
                            ):
                                is_propargyl_alcohol = True

                    # Check if this is a propargylation reaction
                    if (
                        ketone_in_reactants
                        and is_propargyl_alcohol
                        and (is_relevant_reaction or has_propargyl_reagent)
                    ):
                        print(f"Found propargylation of ketone: {rsmi}")
                        has_propargylation = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Propargylation of ketone strategy detected: {has_propargylation}")
    return has_propargylation
