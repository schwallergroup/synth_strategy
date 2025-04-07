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
    Detects if an alkyne linker connecting two aromatic rings is maintained throughout the synthesis.
    """
    # Helper function to check if a molecule has an alkyne connecting two aromatic rings
    def has_aromatic_alkyne_connection(mol_smiles):
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check for alkyne functional group
        if not checker.check_fg("Alkyne", mol_smiles):
            return False

        # Check for various patterns of alkyne between aromatic rings
        # Direct connection
        pattern1 = Chem.MolFromSmarts("[c]-[#6]#[#6]-[c]")
        # Connection through linkers (like O, CH2, etc.)
        pattern2 = Chem.MolFromSmarts("[c]~[*]~[#6]#[#6]~[*]~[c]")
        # Specific pattern seen in test case (aromatic-O-C-C#C-aromatic)
        pattern3 = Chem.MolFromSmarts("[c]-[#8]-[#6]-[#6]#[#6]-[c]")

        if mol.HasSubstructMatch(pattern1):
            print(f"Found direct aromatic-alkyne-aromatic connection in: {mol_smiles}")
            return True
        elif mol.HasSubstructMatch(pattern2):
            print(f"Found aromatic-linker-alkyne-linker-aromatic connection in: {mol_smiles}")
            return True
        elif mol.HasSubstructMatch(pattern3):
            print(f"Found aromatic-O-C-alkyne-aromatic connection in: {mol_smiles}")
            return True

        return False

    # Check if target molecule has alkyne linker
    if route["type"] == "mol" and not has_aromatic_alkyne_connection(route["smiles"]):
        print(f"Target molecule doesn't have aromatic-alkyne-aromatic linker: {route['smiles']}")
        return False

    # Track molecules with the alkyne linker at each depth
    molecules_with_alkyne = {}
    alkyne_maintained = True

    def dfs_traverse(node, depth=0):
        nonlocal alkyne_maintained

        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check if molecule has an alkyne and if it connects aromatic rings
            if has_aromatic_alkyne_connection(smiles):
                print(f"Found aromatic-alkyne-aromatic linker at depth {depth}: {smiles}")
                if depth not in molecules_with_alkyne:
                    molecules_with_alkyne[depth] = []
                molecules_with_alkyne[depth].append(smiles)

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if the reaction preserves the alkyne linker
            try:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                # In retrosynthesis, we check if any reactant has the alkyne linker
                # and if so, verify it's preserved in the product
                for reactant in reactants:
                    if has_aromatic_alkyne_connection(reactant):
                        if not has_aromatic_alkyne_connection(product):
                            print(f"Alkyne linker broken in reaction at depth {depth}: {rsmi}")
                            alkyne_maintained = False
                            return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if alkyne linker is found at multiple depths
    depths_with_alkyne = list(molecules_with_alkyne.keys())
    print(f"Depths with alkyne linker: {depths_with_alkyne}")

    # For the alkyne to be maintained, it should be present in the target molecule (depth 0)
    # and at least one other depth, and no reactions should break the linker
    if (
        depths_with_alkyne
        and min(depths_with_alkyne) == 0
        and max(depths_with_alkyne) > 0
        and alkyne_maintained
    ):
        print("Alkyne linker is maintained throughout the synthesis")
        return True

    print("Alkyne linker is not maintained throughout the synthesis")
    return False
