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
    This function detects phthalimide protection of primary amine in the synthetic route.
    """
    protection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal protection_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction SMILES: {rsmi}")

            # Check if this is a phthalimide protection or deprotection reaction
            if checker.check_reaction("Phthalic anhydride to phthalimide", rsmi):
                print("Phthalimide formation reaction detected")
                protection_detected = True
                return

            if checker.check_reaction("Phthalimide deprotection", rsmi):
                print("Phthalimide deprotection reaction detected")
                protection_detected = True
                return

            # If not directly detected, check for the functional group transformation
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for primary amine in reactants
                has_primary_amine = False
                for reactant in reactants:
                    if checker.check_fg("Primary amine", reactant):
                        print(f"Primary amine found in reactant: {reactant}")
                        has_primary_amine = True
                        break

                # Check for phthalimide in product
                if has_primary_amine:
                    prod_mol = Chem.MolFromSmiles(product)

                    # Check for N-substituted phthalimide structure
                    if checker.check_fg("Unsubstituted dicarboximide", product):
                        print("Dicarboximide group found in product")
                        protection_detected = True
                        return

                    # Check for phthalimide ring structure
                    if prod_mol and prod_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("O=C1NC(=O)c2ccccc21")
                    ):
                        print("Phthalimide ring structure found in product")
                        protection_detected = True
                        return

                # Check for phthalimide in reactants and primary amine in product (deprotection)
                has_phthalimide = False
                for reactant in reactants:
                    if checker.check_fg("Unsubstituted dicarboximide", reactant) or (
                        Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(
                            Chem.MolFromSmarts("O=C1NC(=O)c2ccccc21")
                        )
                    ):
                        print(f"Phthalimide found in reactant: {reactant}")
                        has_phthalimide = True
                        break

                if has_phthalimide and checker.check_fg("Primary amine", product):
                    print(
                        "Primary amine found in product after phthalimide deprotection"
                    )
                    protection_detected = True
                    return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Protection detected: {protection_detected}")
    return protection_detected
