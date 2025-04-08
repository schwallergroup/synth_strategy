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
    Detects if the synthesis includes an O-demethylation step (conversion of methoxy to hydroxyl).
    """
    demethylation_detected = False

    def dfs_traverse(node):
        nonlocal demethylation_detected
        if demethylation_detected:
            return  # Early return if already detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for known demethylation reaction types
            if (
                checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
            ):
                print(f"Detected O-demethylation via reaction type: {rsmi}")
                demethylation_detected = True
                return

            # If not identified by reaction type, check for the functional group transformation
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Split reactants and check each one
                reactants = reactants_part.split(".")
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol:
                    # Check for any hydroxyl group in product (phenol or alcohol)
                    has_hydroxyl = (
                        checker.check_fg("Phenol", product_part)
                        or checker.check_fg("Primary alcohol", product_part)
                        or checker.check_fg("Secondary alcohol", product_part)
                        or checker.check_fg("Tertiary alcohol", product_part)
                    )

                    if has_hydroxyl:
                        for reactant in reactants:
                            # Check if reactant has a methoxy group
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Check for methoxy group (both aromatic and aliphatic)
                                methoxy_pattern = Chem.MolFromSmarts("[#6]-O[CH3]")
                                if reactant_mol.HasSubstructMatch(methoxy_pattern):
                                    # For atom-mapped reactions, verify the transformation occurs at the same position
                                    if (
                                        "[" in rsmi and ":" in rsmi
                                    ):  # Check if reaction is atom-mapped
                                        # Get atom indices of methoxy in reactant
                                        matches = reactant_mol.GetSubstructMatches(methoxy_pattern)
                                        for match in matches:
                                            carbon_idx = match[0]  # Carbon attached to oxygen
                                            # Check if this carbon has a hydroxyl in the product
                                            carbon_map_num = None
                                            for atom in reactant_mol.GetAtoms():
                                                if atom.GetIdx() == carbon_idx:
                                                    if atom.HasProp("molAtomMapNumber"):
                                                        carbon_map_num = atom.GetProp(
                                                            "molAtomMapNumber"
                                                        )
                                                    break

                                            if carbon_map_num:
                                                # Look for this mapped atom in product with hydroxyl
                                                hydroxyl_pattern = Chem.MolFromSmarts("[#6]-[OH]")
                                                for atom in product_mol.GetAtoms():
                                                    if (
                                                        atom.HasProp("molAtomMapNumber")
                                                        and atom.GetProp("molAtomMapNumber")
                                                        == carbon_map_num
                                                    ):
                                                        # Check if this atom is part of a hydroxyl group
                                                        for (
                                                            match
                                                        ) in product_mol.GetSubstructMatches(
                                                            hydroxyl_pattern
                                                        ):
                                                            if (
                                                                atom.GetIdx() == match[0]
                                                            ):  # Carbon attached to hydroxyl
                                                                print(
                                                                    f"Detected O-demethylation via atom mapping: {rsmi}"
                                                                )
                                                                demethylation_detected = True
                                                                return
                                    else:
                                        # For non-atom-mapped reactions, just confirm the presence of both groups
                                        print(
                                            f"Detected O-demethylation via functional group analysis: {rsmi}"
                                        )
                                        demethylation_detected = True
                                        return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return demethylation_detected
