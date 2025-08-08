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
    Detects if a brominated pyridine scaffold is preserved throughout the synthesis.
    """
    # Track depths with brominated pyridine
    depths_with_bromopyridine = set()
    mol_depths = set()

    def has_bromopyridine(mol_smiles):
        """Check if a molecule contains a brominated pyridine scaffold"""
        # First check if it has both a pyridine ring and an aromatic halide
        if not checker.check_ring("pyridine", mol_smiles):
            print(f"No pyridine ring found in {mol_smiles}")
            return False

        # Check for bromine specifically
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            print(f"Could not parse molecule: {mol_smiles}")
            return False

        # Check for aromatic bromine
        has_aromatic_br = False
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "Br" and atom.GetIsAromatic():
                has_aromatic_br = True
                break
            # Check if it's a bromine attached to an aromatic atom
            if atom.GetSymbol() == "Br":
                for bond in atom.GetBonds():
                    other_atom = mol.GetAtomWithIdx(bond.GetOtherAtomIdx(atom.GetIdx()))
                    if other_atom.GetIsAromatic():
                        has_aromatic_br = True
                        break

        if not has_aromatic_br and not checker.check_fg("Aromatic halide", mol_smiles):
            print(f"No aromatic bromine found in {mol_smiles}")
            return False

        # Use multiple SMARTS patterns to catch different bromopyridine arrangements
        bromopyridine_patterns = [
            "c1ncccc1[Br]",  # 2-bromopyridine
            "c1nccc([Br])c1",  # 3-bromopyridine
            "c1ncc([Br])cc1",  # 4-bromopyridine
            "c1nc([Br])ccc1",  # 5-bromopyridine
            "c1nc(ccc1[Br])",  # 6-bromopyridine
        ]

        for pattern in bromopyridine_patterns:
            patt = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(patt):
                print(f"Found brominated pyridine via SMARTS in {mol_smiles}")
                return True

        # If SMARTS patterns didn't match, try a more detailed approach
        try:
            # Get pyridine ring atoms
            pyridine_indices = checker.get_ring_atom_indices("pyridine", mol_smiles)
            if not pyridine_indices:
                print(f"No pyridine indices found in {mol_smiles}")
                return False

            pyridine_atoms = set()
            # Handle different possible return formats
            if isinstance(pyridine_indices, list):
                for ring_atoms in pyridine_indices:
                    if isinstance(ring_atoms, tuple) and len(ring_atoms) > 0:
                        if isinstance(ring_atoms[0], tuple):
                            for atom_idx in ring_atoms[0]:
                                pyridine_atoms.add(atom_idx)
                        elif isinstance(ring_atoms[0], int):
                            for atom_idx in ring_atoms:
                                pyridine_atoms.add(atom_idx)
            elif isinstance(pyridine_indices, int):
                # Handle case where a single index is returned
                pyridine_atoms.add(pyridine_indices)

            # Check if any bromine is attached to the pyridine ring
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "Br":
                    for bond in atom.GetBonds():
                        other_atom_idx = bond.GetOtherAtomIdx(atom.GetIdx())
                        if other_atom_idx in pyridine_atoms:
                            print(f"Found brominated pyridine in {mol_smiles}")
                            return True

            print(f"No brominated pyridine found in {mol_smiles}")
            return False

        except Exception as e:
            print(f"Error checking for brominated pyridine: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_bromopyridine, mol_depths

        if node["type"] == "mol":
            mol_depths.add(depth)
            if has_bromopyridine(node["smiles"]):
                depths_with_bromopyridine.add(depth)

                # Check if this is a starting material
                if node.get("in_stock", False):
                    print(f"Found brominated pyridine in starting material at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if bromopyridine is present at all molecule depths
    is_preserved = len(depths_with_bromopyridine) == len(mol_depths) and len(mol_depths) > 0

    print(f"Bromopyridine scaffold preservation: {is_preserved}")
    print(f"Depths with bromopyridine: {depths_with_bromopyridine}")
    print(f"Molecule depths: {mol_depths}")

    return is_preserved
