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
    This function detects the use of orthogonal protection strategy
    with both Boc and Cbz protecting groups in the same synthesis.
    """
    # Track protection and deprotection for both groups
    has_boc_group = False
    has_cbz_group = False
    has_boc_protection = False
    has_boc_deprotection = False
    has_cbz_protection = False
    has_cbz_deprotection = False

    # Create SMARTS patterns for Boc and Cbz groups
    boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
    cbz_pattern = Chem.MolFromSmarts("c1ccccc1COC(=O)[N]")

    def dfs_traverse(node):
        nonlocal has_boc_group, has_cbz_group, has_boc_protection, has_boc_deprotection, has_cbz_protection, has_cbz_deprotection

        # Check molecule nodes for protecting groups
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check for Boc and Cbz groups in molecules
            if checker.check_fg("Carbamic ester", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)

                # Check for Boc group (tert-butyloxycarbonyl)
                if mol and mol.HasSubstructMatch(boc_pattern):
                    print(f"Found Boc protecting group in molecule: {mol_smiles}")
                    has_boc_group = True

                # Check for Cbz group (benzyloxycarbonyl)
                if mol and mol.HasSubstructMatch(cbz_pattern):
                    print(f"Found Cbz protecting group in molecule: {mol_smiles}")
                    has_cbz_group = True

        # Check reaction nodes for protection/deprotection reactions
        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for Boc protection reactions
            if (
                checker.check_reaction("Boc amine protection", rxn_smiles)
                or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rxn_smiles)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rxn_smiles)
                or checker.check_reaction("Boc amine protection of secondary amine", rxn_smiles)
                or checker.check_reaction("Boc amine protection of primary amine", rxn_smiles)
            ):
                print(f"Found Boc protection reaction: {rxn_smiles}")
                has_boc_protection = True

            # Check for Boc deprotection reactions
            if (
                checker.check_reaction("Boc amine deprotection", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection of guanidine", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rxn_smiles)
                or checker.check_reaction("Tert-butyl deprotection of amine", rxn_smiles)
            ):
                print(f"Found Boc deprotection reaction: {rxn_smiles}")
                has_boc_deprotection = True

            # Check for Cbz protection/deprotection
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # For Cbz protection, look for reactions that form carbamic esters with benzyl groups
            if checker.check_fg("Carbamic ester", product):
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(cbz_pattern):
                    print(f"Found Cbz protection reaction: {rxn_smiles}")
                    has_cbz_protection = True

            # For Cbz deprotection, check for hydrogenolysis of carbamic esters
            if (
                checker.check_reaction("Carboxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction("Hydroxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction("Hydrogenolysis of amides/imides/carbamates", rxn_smiles)
            ):
                # Verify it's specifically removing a Cbz group
                for reactant in reactants:
                    if reactant and checker.check_fg("Carbamic ester", reactant):
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol and react_mol.HasSubstructMatch(cbz_pattern):
                            print(f"Found Cbz deprotection reaction: {rxn_smiles}")
                            has_cbz_deprotection = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Boc group: {has_boc_group}, Cbz group: {has_cbz_group}")
    print(f"Boc protection: {has_boc_protection}, Boc deprotection: {has_boc_deprotection}")
    print(f"Cbz protection: {has_cbz_protection}, Cbz deprotection: {has_cbz_deprotection}")

    # The orthogonal protection strategy is present if both protecting groups are used
    # We don't strictly require seeing the protection/deprotection reactions
    return has_boc_group and has_cbz_group
