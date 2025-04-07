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
    This function detects if the synthesis uses a Boc-protected piperazine scaffold.
    """
    boc_piperazine_present = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_piperazine_present

        if node["type"] == "mol":
            try:
                mol_smiles = node["smiles"]

                # Check if molecule contains piperazine ring
                has_piperazine = checker.check_ring("piperazine", mol_smiles)

                # Check if molecule contains Boc group
                has_boc = checker.check_fg("Boc", mol_smiles)

                # Simple check for Boc-protected piperazine
                if has_piperazine and has_boc:
                    boc_piperazine_present = True
                    print(
                        f"Boc-protected piperazine molecule detected at depth {depth}: {mol_smiles}"
                    )

                    # Create molecule to check if Boc is attached to piperazine
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Get atom indices for piperazine ring
                        piperazine_indices = checker.get_ring_atom_indices("piperazine", mol_smiles)
                        if piperazine_indices:
                            # Get nitrogen atoms in piperazine
                            piperazine_n_indices = []
                            for ring_atoms in piperazine_indices:
                                for atom_idx in ring_atoms:
                                    atom = mol.GetAtomWithIdx(atom_idx)
                                    if atom.GetAtomicNum() == 7:  # Nitrogen
                                        piperazine_n_indices.append(atom_idx)

                            # Get Boc group atoms
                            boc_indices = checker.get_fg_atom_indices("Boc", mol_smiles)

                            # Check if any Boc group is attached to a piperazine nitrogen
                            if boc_indices and piperazine_n_indices:
                                for boc_atoms in boc_indices:
                                    for atom_idx in boc_atoms:
                                        atom = mol.GetAtomWithIdx(atom_idx)
                                        for neighbor in atom.GetNeighbors():
                                            if neighbor.GetIdx() in piperazine_n_indices:
                                                print(
                                                    f"Confirmed: Boc group is directly attached to piperazine nitrogen at depth {depth}"
                                                )
                                                break
            except Exception as e:
                print(f"Error processing molecule {mol_smiles}: {e}")

        elif node["type"] == "reaction":
            try:
                if "metadata" in node and "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]

                    # Check for Boc protection/deprotection reactions
                    is_boc_protection = checker.check_reaction("Boc amine protection", rsmi)
                    is_boc_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)

                    if is_boc_protection or is_boc_deprotection:
                        # Extract reactants and product
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant or product contains piperazine
                        for mol_smiles in reactants + [product]:
                            if checker.check_ring("piperazine", mol_smiles):
                                boc_piperazine_present = True
                                print(
                                    f"Boc-{'protection' if is_boc_protection else 'deprotection'} of piperazine detected at depth {depth}"
                                )
                                break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return boc_piperazine_present
