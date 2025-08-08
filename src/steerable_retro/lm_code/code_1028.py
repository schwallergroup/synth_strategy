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
    This function detects if a halogenated aromatic scaffold (specifically
    3,5-dichloro-4-fluoro pattern) is maintained throughout the synthesis.
    """
    all_nodes_have_pattern = True

    def dfs_traverse(node):
        nonlocal all_nodes_have_pattern

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Check for aromatic ring with 3,5-dichloro-4-fluoro pattern
                # This pattern allows for other substituents to be present
                pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)c(F)c1")

                # Alternative pattern to catch different atom ordering
                alt_pattern = Chem.MolFromSmarts("c1c(Cl)c(F)c(Cl)cc1")

                # Check if molecule has aromatic halides
                has_aromatic_halide = checker.check_fg("Aromatic halide", mol_smiles)

                # Check for the specific pattern
                if not (mol.HasSubstructMatch(pattern) or mol.HasSubstructMatch(alt_pattern)):
                    # If it's a small molecule like "N" or a reagent, don't consider it
                    if len(mol.GetAtoms()) > 5 and has_aromatic_halide:
                        all_nodes_have_pattern = False
                        print(f"Molecule does not have the halogenated pattern: {mol_smiles}")
                else:
                    print(f"Molecule has the halogenated pattern: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Halogenated aromatic scaffold maintained throughout: {all_nodes_have_pattern}")
    return all_nodes_have_pattern
