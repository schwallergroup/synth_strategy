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
    Detects if the synthesis uses late-stage functional group modifications
    without changing the carbon skeleton in the final steps.
    """
    late_stage_mods = 0
    total_steps = 0

    # List of functional groups to check
    fg_list = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Phenol",
        "Carboxylic acid",
        "Ester",
        "Aldehyde",
        "Ketone",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Amide",
        "Nitrile",
        "Halide",
        "Nitro group",
        "Ether",
        "Acyl halide",
        "Anhydride",
        "Sulfonamide",
    ]

    # List of reaction types that typically involve functional group modifications
    fg_mod_reactions = [
        "Esterification of Carboxylic Acids",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Reduction of aldehydes and ketones to alcohols",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Alcohol protection with silyl ethers",
        "Alcohol deprotection from silyl ethers",
        "Boc amine protection",
        "Boc amine deprotection",
        "Oxidation of aldehydes to carboxylic acids",
        "Reduction of ester to primary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Schotten-Baumann to ester",
        "Acetic anhydride and alcohol to ester",
        "Acylation of secondary amines",
        "Acylation of primary amines",
        "Acylation of olefines by aldehydes",
        "Acylation of secondary amines with anhydrides",
    ]

    def count_carbon_atoms(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
        return 0

    def get_mapped_atoms(smiles):
        """Extract atom mapping dictionary from SMILES"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}

        atom_map = {}
        for atom in mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                atom_map[map_num] = atom.GetIdx()
        return atom_map

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_mods, total_steps

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            total_steps += 1

            # For late stage (depth < 4), check for functional group modifications
            if depth < 4:
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a known functional group modification reaction
                is_fg_mod_reaction = False
                for rxn_type in fg_mod_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction at depth {depth}")
                        is_fg_mod_reaction = True
                        break

                # Check carbon skeleton preservation
                product_carbons = count_carbon_atoms(product)
                reactant_carbons_list = [count_carbon_atoms(r) for r in reactants]

                # Find the main reactant (the one with the most carbon atoms)
                if reactant_carbons_list:
                    main_reactant_idx = reactant_carbons_list.index(max(reactant_carbons_list))
                    main_reactant = reactants[main_reactant_idx]
                    main_reactant_carbons = reactant_carbons_list[main_reactant_idx]

                    # Carbon skeleton is preserved if carbon count difference is small
                    # For larger molecules, allow slightly larger differences
                    max_carbon_diff = 2
                    if main_reactant_carbons > 20:
                        max_carbon_diff = 3

                    carbon_preserved = (
                        abs(main_reactant_carbons - product_carbons) <= max_carbon_diff
                    )

                    if carbon_preserved:
                        # Check for functional group changes using the checker

                        # Check functional groups in reactant and product
                        reactant_fgs = []
                        product_fgs = []

                        for fg in fg_list:
                            if checker.check_fg(fg, main_reactant):
                                reactant_fgs.append(fg)
                            if checker.check_fg(fg, product):
                                product_fgs.append(fg)

                        # Check if functional groups have changed
                        fg_diff = set(reactant_fgs).symmetric_difference(set(product_fgs))

                        if fg_diff or is_fg_mod_reaction:
                            print(
                                f"Late-stage functional group modification detected at depth {depth}"
                            )
                            print(f"Reactant FGs: {reactant_fgs}")
                            print(f"Product FGs: {product_fgs}")
                            print(f"FG differences: {fg_diff}")
                            late_stage_mods += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Total late-stage modifications: {late_stage_mods}")
    # If at least 1 late-stage modification is found
    return late_stage_mods >= 1
