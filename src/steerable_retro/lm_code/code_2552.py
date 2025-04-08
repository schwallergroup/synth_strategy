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
    This function detects conversion of alcohol to alkyl chloride in the synthetic route.
    """
    transformation_detected = False

    def dfs_traverse(node):
        nonlocal transformation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an alcohol to chloride reaction using the checker
                alcohol_to_chloride_reactions = [
                    "Alcohol to chloride_sulfonyl chloride",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_CHCl3",
                    "Alcohol to chloride_CH2Cl2",
                    "Alcohol to chloride_PCl5_ortho",
                    "Alcohol to chloride_POCl3_ortho",
                    "Alcohol to chloride_POCl3_para",
                    "Alcohol to chloride_POCl3",
                    "Alcohol to chloride_HCl",
                    "Alcohol to chloride_Salt",
                    "Alcohol to chloride_Other",
                ]

                for reaction_type in alcohol_to_chloride_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} reaction")
                        transformation_detected = True
                        break

                # If no specific reaction type was detected, check for alcohol in reactants and chloride in product
                if not transformation_detected:
                    try:
                        # Check for various alcohol types in reactants
                        alcohol_types = [
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Aromatic alcohol",
                        ]
                        has_alcohol_reactant = False

                        for reactant in reactants:
                            for alcohol_type in alcohol_types:
                                if checker.check_fg(alcohol_type, reactant):
                                    print(f"Found {alcohol_type} in reactant: {reactant}")
                                    has_alcohol_reactant = True
                                    break
                            if has_alcohol_reactant:
                                break

                        # Check for chloride in product
                        has_chloride_product = False
                        if (
                            checker.check_fg("Primary halide", product)
                            or checker.check_fg("Secondary halide", product)
                            or checker.check_fg("Tertiary halide", product)
                        ):
                            # Verify it's specifically a chloride
                            prod_mol = Chem.MolFromSmiles(product)
                            if prod_mol and prod_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[C]-[Cl]")
                            ):
                                print(f"Found chloride in product: {product}")
                                has_chloride_product = True

                        if has_alcohol_reactant and has_chloride_product:
                            print(
                                "Alcohol to chloride conversion detected through functional group analysis"
                            )
                            transformation_detected = True
                    except Exception as e:
                        print(f"Error in functional group analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_detected
