#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects if the final product contains multiple sulfonyl groups.
    """
    # The final product is the root node in retrosynthetic analysis
    final_product_smiles = None

    if route["type"] == "mol":
        final_product_smiles = route["smiles"]
        print(f"Found final product: {final_product_smiles}")

    result = False
    if final_product_smiles:
        # Count sulfonyl groups using the checker functions
        sulfonyl_count = 0

        # Debug: Print all functional groups in the molecule
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol:
            print(f"Molecule has {mol.GetNumAtoms()} atoms")

        # Check for sulfones
        if checker.check_fg("Sulfone", final_product_smiles):
            sulfone_indices = checker.get_fg_atom_indices("Sulfone", final_product_smiles)
            if sulfone_indices:
                sulfonyl_count += len(sulfone_indices)
                print(f"Found {len(sulfone_indices)} sulfone groups")

        # Check for sulfonamides
        if checker.check_fg("Sulfonamide", final_product_smiles):
            sulfonamide_indices = checker.get_fg_atom_indices("Sulfonamide", final_product_smiles)
            if sulfonamide_indices:
                sulfonyl_count += len(sulfonamide_indices)
                print(f"Found {len(sulfonamide_indices)} sulfonamide groups")

        # Check for sulfonyl halides
        if checker.check_fg("Sulfonyl halide", final_product_smiles):
            sulfonyl_halide_indices = checker.get_fg_atom_indices(
                "Sulfonyl halide", final_product_smiles
            )
            if sulfonyl_halide_indices:
                sulfonyl_count += len(sulfonyl_halide_indices)
                print(f"Found {len(sulfonyl_halide_indices)} sulfonyl halide groups")

        # Check for sulfonic acids
        if checker.check_fg("Sulfonic acid", final_product_smiles):
            sulfonic_acid_indices = checker.get_fg_atom_indices(
                "Sulfonic acid", final_product_smiles
            )
            if sulfonic_acid_indices:
                sulfonyl_count += len(sulfonic_acid_indices)
                print(f"Found {len(sulfonic_acid_indices)} sulfonic acid groups")

        # Check for sulfonates
        if checker.check_fg("Sulfonate", final_product_smiles):
            sulfonate_indices = checker.get_fg_atom_indices("Sulfonate", final_product_smiles)
            if sulfonate_indices:
                sulfonyl_count += len(sulfonate_indices)
                print(f"Found {len(sulfonate_indices)} sulfonate groups")

        # Check for sulfoxides
        if checker.check_fg("Sulfoxide", final_product_smiles):
            sulfoxide_indices = checker.get_fg_atom_indices("Sulfoxide", final_product_smiles)
            if sulfoxide_indices:
                sulfonyl_count += len(sulfoxide_indices)
                print(f"Found {len(sulfoxide_indices)} sulfoxide groups")

        # Direct pattern matching as a fallback
        # This is needed because the test case shows the checker might be missing some patterns
        mol = Chem.MolFromSmiles(final_product_smiles)
        if mol:
            # Count S(=O)(=O) patterns directly
            patt = Chem.MolFromSmarts("S(=O)(=O)")
            matches = mol.GetSubstructMatches(patt)
            direct_count = len(matches)
            print(f"Direct S(=O)(=O) pattern matching found {direct_count} instances")

            # If direct matching finds more groups than our checker functions, use that count
            if direct_count > sulfonyl_count:
                print(f"Using direct pattern count: {direct_count} instead of {sulfonyl_count}")
                sulfonyl_count = direct_count

        result = sulfonyl_count >= 2
        print(f"Found {sulfonyl_count} total sulfonyl groups in final product")

    print(f"Multiple sulfonyl groups in final product strategy detected: {result}")
    return result
