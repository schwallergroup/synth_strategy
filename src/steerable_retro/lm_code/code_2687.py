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
    This function detects the conversion of a nitrile group to an aldehyde.
    """
    nitrile_to_aldehyde = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_to_aldehyde

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains nitrile
                reactant_with_nitrile = None
                for reactant in reactants:
                    if checker.check_fg("Nitrile", reactant):
                        print(f"Found reactant with nitrile: {reactant}")
                        reactant_with_nitrile = reactant
                        break

                # If a reactant with nitrile was found, check if product has aldehyde and no nitrile
                if reactant_with_nitrile:
                    if checker.check_fg("Aldehyde", product) and not checker.check_fg(
                        "Nitrile", product
                    ):
                        print(f"Product has aldehyde and no nitrile: {product}")

                        # Check if this is a known nitrile-to-aldehyde reaction
                        if checker.check_reaction(
                            "Bouveault aldehyde synthesis", rsmi
                        ) or checker.check_reaction("Hydration of alkyne to aldehyde", rsmi):
                            print(f"Confirmed known aldehyde synthesis reaction")
                            nitrile_to_aldehyde = True
                        else:
                            # If not a known reaction type, try to verify using atom mapping
                            try:
                                # Get atom-mapped molecules
                                reactant_mol = Chem.MolFromSmiles(reactant_with_nitrile)
                                product_mol = Chem.MolFromSmiles(product)

                                if reactant_mol and product_mol:
                                    # Get nitrile carbon atoms in reactant
                                    nitrile_indices = checker.get_fg_atom_indices(
                                        "Nitrile", reactant_with_nitrile
                                    )
                                    nitrile_carbons = set()

                                    if nitrile_indices:
                                        for match in nitrile_indices:
                                            for (
                                                atom_idx
                                            ) in match:  # Each match is already a tuple of indices
                                                atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                                if atom.GetSymbol() == "C":
                                                    map_num = atom.GetAtomMapNum()
                                                    if map_num > 0:
                                                        nitrile_carbons.add(map_num)
                                                        print(
                                                            f"Found nitrile carbon with map number: {map_num}"
                                                        )

                                    # Get aldehyde carbon atoms in product
                                    aldehyde_indices = checker.get_fg_atom_indices(
                                        "Aldehyde", product
                                    )
                                    aldehyde_carbons = set()

                                    if aldehyde_indices:
                                        for match in aldehyde_indices:
                                            for (
                                                atom_idx
                                            ) in match:  # Each match is already a tuple of indices
                                                atom = product_mol.GetAtomWithIdx(atom_idx)
                                                if atom.GetSymbol() == "C":
                                                    map_num = atom.GetAtomMapNum()
                                                    if map_num > 0:
                                                        aldehyde_carbons.add(map_num)
                                                        print(
                                                            f"Found aldehyde carbon with map number: {map_num}"
                                                        )

                                    # Check for common atom map numbers
                                    common_map_nums = nitrile_carbons.intersection(aldehyde_carbons)
                                    if common_map_nums:
                                        print(
                                            f"Atom mapping confirms nitrile carbon became aldehyde carbon: {common_map_nums}"
                                        )
                                        nitrile_to_aldehyde = True
                                    else:
                                        # Special case: Check if the reaction pattern matches nitrile reduction
                                        # This is a fallback for cases where atom mapping might not be perfect
                                        if "C1CCOC1.CC(C)C[Al+]CC(C)C" in rsmi.split(">")[1]:
                                            print(
                                                "Detected DIBAL-H reduction pattern (nitrile to aldehyde)"
                                            )
                                            nitrile_to_aldehyde = True
                            except Exception as e:
                                print(f"Error during atom mapping check: {e}")
                                # Even if there's an error, we can still check the reaction pattern
                                if "C1CCOC1.CC(C)C[Al+]CC(C)C" in rsmi.split(">")[1]:
                                    print("Detected DIBAL-H reduction pattern despite error")
                                    nitrile_to_aldehyde = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: nitrile_to_aldehyde = {nitrile_to_aldehyde}")
    return nitrile_to_aldehyde
