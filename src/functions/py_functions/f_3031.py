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
    This function detects if the synthetic route uses Boc-protected piperazine.
    """
    has_boc_piperazine = False

    def dfs_traverse(node):
        nonlocal has_boc_piperazine

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if molecule contains both piperazine ring and Boc group
            if checker.check_ring("piperazine", mol_smiles) and checker.check_fg(
                "Boc", mol_smiles
            ):
                # Get the atom indices for both the piperazine ring and Boc group
                piperazine_indices = checker.get_ring_atom_indices(
                    "piperazine", mol_smiles
                )
                boc_indices = checker.get_fg_atom_indices("Boc", mol_smiles)

                # Convert to RDKit molecule to check connectivity
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check if any piperazine nitrogen is connected to the Boc group
                    for ring_atoms in piperazine_indices:
                        # Flatten the tuple of tuples if needed
                        if isinstance(ring_atoms[0], tuple):
                            ring_atoms = [
                                atom for sublist in ring_atoms for atom in sublist
                            ]
                        else:
                            ring_atoms = ring_atoms

                        # Find nitrogen atoms in the piperazine ring
                        piperazine_nitrogens = [
                            idx
                            for idx in ring_atoms
                            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7
                        ]

                        for boc_atoms in boc_indices:
                            # Flatten the tuple of tuples if needed
                            if isinstance(boc_atoms[0], tuple):
                                boc_atoms = [
                                    atom for sublist in boc_atoms for atom in sublist
                                ]
                            else:
                                boc_atoms = boc_atoms

                            # Check if any Boc atom is bonded to any piperazine nitrogen
                            for n_idx in piperazine_nitrogens:
                                for boc_idx in boc_atoms:
                                    bond = mol.GetBondBetweenAtoms(n_idx, boc_idx)
                                    if bond is not None:
                                        print(
                                            f"Boc-protected piperazine detected: {mol_smiles}"
                                        )
                                        has_boc_piperazine = True
                                        return

        elif node["type"] == "reaction":
            # Check if this is a Boc protection or deprotection reaction involving piperazine
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection with Boc anhydride", rsmi
                    )
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction(
                        "Boc amine protection of secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Boc amine protection of primary amine", rsmi
                    )
                ):

                    # Check if reactants or products contain piperazine
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    for reactant in reactants:
                        if checker.check_ring("piperazine", reactant):
                            print(
                                f"Boc protection of piperazine detected in reaction: {rsmi}"
                            )
                            has_boc_piperazine = True
                            return

                    if checker.check_ring("piperazine", product) and checker.check_fg(
                        "Boc", product
                    ):
                        print(f"Boc-protected piperazine produced in reaction: {rsmi}")
                        has_boc_piperazine = True
                        return

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction(
                        "Boc amine deprotection of guanidine", rsmi
                    )
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):

                    reactants = rsmi.split(">")[0].split(".")

                    for reactant in reactants:
                        if checker.check_ring(
                            "piperazine", reactant
                        ) and checker.check_fg("Boc", reactant):
                            print(
                                f"Deprotection of Boc-piperazine detected in reaction: {rsmi}"
                            )
                            has_boc_piperazine = True
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Boc-protected piperazine detected: {has_boc_piperazine}")
    return has_boc_piperazine
