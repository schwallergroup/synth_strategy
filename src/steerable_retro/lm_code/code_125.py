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
    This function detects a sequence of nitrogen functional group interconversions:
    nitro → amine → azide

    In retrosynthesis, we're looking for: azide → amine → nitro
    """
    nitro_to_amine = False
    amine_to_azide = False
    nitro_depth = -1
    amine_depth = -1
    azide_depth = -1

    # Track atom mappings to ensure same nitrogen atom is converted
    nitro_n_mapping = None
    amine_n_mapping = None
    azide_n_mapping = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_to_amine, amine_to_azide, nitro_depth, amine_depth, azide_depth
        nonlocal nitro_n_mapping, amine_n_mapping, azide_n_mapping

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for functional groups in molecules
            if checker.check_fg("Nitro group", mol_smiles):
                print(f"Molecule with nitro group found at depth {depth}: {mol_smiles}")
                if nitro_depth == -1:
                    nitro_depth = depth

            if checker.check_fg("Primary amine", mol_smiles):
                print(f"Molecule with primary amine found at depth {depth}: {mol_smiles}")
                if amine_depth == -1:
                    amine_depth = depth

            if checker.check_fg("Azide", mol_smiles):
                print(f"Molecule with azide found at depth {depth}: {mol_smiles}")
                if azide_depth == -1:
                    azide_depth = depth

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for nitro to amine conversion (in forward direction)
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {depth}")
                    nitro_to_amine = True
                else:
                    # Alternative check for nitro reduction
                    reactant_has_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants if r
                    )
                    product_has_amine = (
                        checker.check_fg("Primary amine", product) if product else False
                    )

                    if reactant_has_nitro and product_has_amine:
                        print(f"Detected nitro to amine conversion at depth {depth}")
                        nitro_to_amine = True

                # Check for amine to azide conversion (in forward direction)
                if (
                    checker.check_reaction("Formation of Azides from halogens", rsmi)
                    or checker.check_reaction("Amine to azide", rsmi)
                    or checker.check_reaction("Alcohol to azide", rsmi)
                ):

                    # Verify amine is converted to azide
                    reactant_has_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants if r
                    )
                    product_has_azide = checker.check_fg("Azide", product) if product else False

                    if reactant_has_amine and product_has_azide:
                        print(f"Found amine to azide conversion at depth {depth}")
                        amine_to_azide = True
                else:
                    # Alternative check for amine to azide conversion
                    reactant_has_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants if r
                    )
                    product_has_azide = checker.check_fg("Azide", product) if product else False

                    if reactant_has_amine and product_has_azide:
                        print(f"Detected amine to azide conversion at depth {depth}")
                        amine_to_azide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Conversion summary:")
    print(f"Nitro to amine: {nitro_to_amine} at depth {nitro_depth}")
    print(f"Amine to azide: {amine_to_azide} at depth {amine_depth}")
    print(f"Azide found at depth: {azide_depth}")

    # In retrosynthesis, the sequence should be: azide → amine → nitro
    # This means azide_depth < amine_depth < nitro_depth
    correct_order = (
        azide_depth != -1
        and amine_depth != -1
        and nitro_depth != -1
        and azide_depth < amine_depth
        and amine_depth < nitro_depth
    )

    print(f"Correct order (azide_depth < amine_depth < nitro_depth): {correct_order}")
    print(f"Depths - Azide: {azide_depth}, Amine: {amine_depth}, Nitro: {nitro_depth}")

    # Return True if both conversions were found and all three functional groups appear in the correct order
    return nitro_to_amine and amine_to_azide and correct_order
