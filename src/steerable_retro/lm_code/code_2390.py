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
    This function detects a strategy involving C-C bond disconnection
    that produces an aromatic aldehyde fragment.
    """
    has_cc_disconnection = False
    has_aromatic_aldehyde_fragment = False

    def has_carbon_carbon_bond(mol1, mol2, product_mol):
        """Check if there's a C-C bond between fragments in the product"""
        if mol1 is None or mol2 is None or product_mol is None:
            return False

        # Get atom maps from reactants
        atom_maps1 = {
            atom.GetIdx(): atom.GetAtomMapNum()
            for atom in mol1.GetAtoms()
            if atom.GetAtomMapNum() > 0
        }
        atom_maps2 = {
            atom.GetIdx(): atom.GetAtomMapNum()
            for atom in mol2.GetAtoms()
            if atom.GetAtomMapNum() > 0
        }

        # Get carbon atoms from both reactants
        carbon_maps1 = {
            idx: map_num
            for idx, map_num in atom_maps1.items()
            if mol1.GetAtomWithIdx(idx).GetSymbol() == "C"
        }
        carbon_maps2 = {
            idx: map_num
            for idx, map_num in atom_maps2.items()
            if mol2.GetAtomWithIdx(idx).GetSymbol() == "C"
        }

        # Check if any carbon in mol1 is bonded to any carbon in mol2 in the product
        for _, map1 in carbon_maps1.items():
            for _, map2 in carbon_maps2.items():
                # Find corresponding atoms in product
                prod_atoms = [
                    (atom.GetIdx(), atom.GetAtomMapNum())
                    for atom in product_mol.GetAtoms()
                    if atom.GetAtomMapNum() == map1 or atom.GetAtomMapNum() == map2
                ]

                if len(prod_atoms) == 2:
                    idx1, _ = prod_atoms[0]
                    idx2, _ = prod_atoms[1]
                    bond = product_mol.GetBondBetweenAtoms(idx1, idx2)
                    if bond is not None:
                        return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_cc_disconnection, has_aromatic_aldehyde_fragment

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, we're looking for reactions where one product
                # is formed from multiple reactants (C-C bond formation in forward direction)
                if len(reactants) > 1 and "." not in product:
                    print(f"Potential C-C disconnection at depth {depth}: {rsmi}")

                    # Check if one of the reactants is an aromatic aldehyde
                    for reactant in reactants:
                        if checker.check_fg("Aldehyde", reactant) and checker.check_ring(
                            "benzene", reactant
                        ):
                            has_aromatic_aldehyde_fragment = True
                            print(
                                f"Detected aromatic aldehyde fragment at depth {depth}: {reactant}"
                            )

                    # Verify C-C bond disconnection by checking pairs of reactants
                    product_mol = Chem.MolFromSmiles(product)
                    for i in range(len(reactants)):
                        for j in range(i + 1, len(reactants)):
                            mol1 = Chem.MolFromSmiles(reactants[i])
                            mol2 = Chem.MolFromSmiles(reactants[j])
                            if has_carbon_carbon_bond(mol1, mol2, product_mol):
                                has_cc_disconnection = True
                                print(f"Confirmed C-C bond disconnection at depth {depth}")
                                break
                        if has_cc_disconnection:
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key features are detected
    return has_cc_disconnection and has_aromatic_aldehyde_fragment
