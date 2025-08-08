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
    This function detects if the final product contains a morpholine amide structure.
    """
    final_product_has_morpholine_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_morpholine_amide

        if node["type"] == "mol" and depth == 0:  # Final product
            if "smiles" in node:
                mol_smiles = node["smiles"]
                print(f"Analyzing final product: {mol_smiles}")

                # Create RDKit molecule for analysis
                mol = Chem.MolFromSmiles(mol_smiles)
                if not mol:
                    print("Failed to create RDKit molecule")
                    return

                # Check if morpholine ring is present
                if checker.check_ring("morpholine", mol_smiles):
                    print(f"Found morpholine ring in molecule")

                    # Check for all types of amides
                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
                    for amide_type in amide_types:
                        if checker.check_fg(amide_type, mol_smiles):
                            print(f"Found {amide_type} in molecule")

                            # Get the atom indices for morpholine ring and amide
                            morpholine_indices = checker.get_ring_atom_indices(
                                "morpholine", mol_smiles
                            )
                            amide_indices = checker.get_fg_atom_indices(amide_type, mol_smiles)

                            print(f"Morpholine indices: {morpholine_indices}")
                            print(f"Amide indices: {amide_indices}")

                            # Check each morpholine ring against each amide group
                            for morph_atoms in morpholine_indices:
                                morph_n_atoms = [
                                    idx
                                    for idx in morph_atoms
                                    if mol.GetAtomWithIdx(idx).GetSymbol() == "N"
                                ]
                                print(f"Morpholine nitrogen atoms: {morph_n_atoms}")

                                for amide_atoms in amide_indices:
                                    print(f"Checking amide atoms: {amide_atoms}")

                                    # Find the nitrogen and carbonyl carbon in the amide
                                    amide_n_atoms = [
                                        idx
                                        for idx in amide_atoms
                                        if mol.GetAtomWithIdx(idx).GetSymbol() == "N"
                                    ]
                                    amide_c_atoms = [
                                        idx
                                        for idx in amide_atoms
                                        if mol.GetAtomWithIdx(idx).GetSymbol() == "C"
                                    ]

                                    print(f"Amide nitrogen atoms: {amide_n_atoms}")
                                    print(f"Amide carbon atoms: {amide_c_atoms}")

                                    # Check if any morpholine nitrogen is also an amide nitrogen
                                    shared_n_atoms = set(morph_n_atoms).intersection(
                                        set(amide_n_atoms)
                                    )
                                    if shared_n_atoms:
                                        print(f"Found shared nitrogen atoms: {shared_n_atoms}")
                                        final_product_has_morpholine_amide = True
                                        return

                                    # Check if morpholine nitrogen is connected to amide carbonyl carbon
                                    for n_idx in morph_n_atoms:
                                        for c_idx in amide_c_atoms:
                                            if mol.GetBondBetweenAtoms(n_idx, c_idx) is not None:
                                                print(
                                                    f"Found morpholine N ({n_idx}) connected to amide carbonyl C ({c_idx})"
                                                )
                                                final_product_has_morpholine_amide = True
                                                return

                                    # Check if any atom in the amide is connected to any atom in the morpholine
                                    for amide_atom in amide_atoms:
                                        for morph_atom in morph_atoms:
                                            if (
                                                mol.GetBondBetweenAtoms(amide_atom, morph_atom)
                                                is not None
                                            ):
                                                print(
                                                    f"Found morpholine atom ({morph_atom}) connected to amide atom ({amide_atom})"
                                                )
                                                # Now check if this forms a morpholine amide pattern
                                                if (
                                                    mol.GetAtomWithIdx(morph_atom).GetSymbol()
                                                    == "N"
                                                    or mol.GetAtomWithIdx(amide_atom).GetSymbol()
                                                    == "N"
                                                ):
                                                    final_product_has_morpholine_amide = True
                                                    return

                # Additional check for morpholine amide pattern
                # Check if the molecule contains a morpholine connected to an amide
                if "NC3CCOCC3" in mol_smiles or "N1CCOCC1" in mol_smiles:
                    print("Found morpholine pattern in SMILES")
                    for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]:
                        if checker.check_fg(amide_type, mol_smiles):
                            print(f"Found {amide_type} with morpholine pattern")
                            final_product_has_morpholine_amide = True
                            return

                # Direct check for C(=O)NC3CCOCC3 pattern (morpholine amide)
                if "C(=O)NC3CCOCC3" in mol_smiles or "C(=O)N1CCOCC1" in mol_smiles:
                    print("Found direct morpholine amide pattern in SMILES")
                    final_product_has_morpholine_amide = True
                    return

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {final_product_has_morpholine_amide}")
    return final_product_has_morpholine_amide
