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
    This function detects if the synthesis preserves a 2,4-difluoro-anisole core
    throughout the synthesis.
    """
    # Track if we find the core in all molecules
    all_molecules_have_core = True
    molecule_count = 0

    def dfs_traverse(node):
        nonlocal all_molecules_have_core, molecule_count

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Skip simple reagents/solvents that wouldn't contain the core
            if mol_smiles in ["O=C=O", "O=C(Cl)C(=O)Cl", "N"]:
                print(f"Skipping reagent/solvent: {mol_smiles}")
                return

            molecule_count += 1
            try:
                mol = Chem.MolFromSmiles(mol_smiles)
                if not mol:
                    print(f"Failed to parse molecule: {mol_smiles}")
                    all_molecules_have_core = False
                    return

                # Check for benzene ring
                has_benzene = checker.check_ring("benzene", mol_smiles)

                if not has_benzene:
                    print(f"Molecule without benzene ring: {mol_smiles}")
                    all_molecules_have_core = False
                    return

                # Check for 2,4-difluoro-anisole core or its phenol precursor
                # Using substructure matching for more reliable detection
                anisole_pattern = Chem.MolFromSmarts("COc1cc(F)cc(F)c1")
                phenol_pattern = Chem.MolFromSmarts("Oc1cc(F)cc(F)c1")

                has_anisole_core = mol.HasSubstructMatch(anisole_pattern)
                has_phenol_core = mol.HasSubstructMatch(phenol_pattern)

                # Also check for variations where other groups are attached to the core
                # This handles cases like the cyano or amide derivatives
                anisole_core_with_substituent = Chem.MolFromSmarts(
                    "COc1cc(F)c([*])c(F)c1"
                )
                phenol_core_with_substituent = Chem.MolFromSmarts(
                    "Oc1cc(F)c([*])c(F)c1"
                )

                has_anisole_with_substituent = mol.HasSubstructMatch(
                    anisole_core_with_substituent
                )
                has_phenol_with_substituent = mol.HasSubstructMatch(
                    phenol_core_with_substituent
                )

                has_core = (
                    has_anisole_core
                    or has_phenol_core
                    or has_anisole_with_substituent
                    or has_phenol_with_substituent
                )

                if not has_core:
                    print(f"Molecule without core or valid precursor: {mol_smiles}")
                    all_molecules_have_core = False
                else:
                    if has_anisole_core or has_anisole_with_substituent:
                        print(f"Molecule with core: {mol_smiles}")
                    else:
                        print(f"Molecule with core precursor (phenol): {mol_smiles}")
            except Exception as e:
                print(f"Error processing molecule {mol_smiles}: {str(e)}")
                all_molecules_have_core = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only return true if we actually checked some molecules and they all had the core
    return all_molecules_have_core and molecule_count > 0
