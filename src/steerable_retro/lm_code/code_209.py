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
    This function detects a synthetic strategy where a heterocyclic scaffold
    (like benzothiophene) is preserved throughout the synthesis while functional
    groups are modified.
    """
    # Track if heterocyclic scaffold is present in all molecules
    all_mols_have_scaffold = True
    mol_count = 0

    # List of heterocyclic scaffolds to check
    heterocyclic_scaffolds = [
        "benzothiophene",
        "benzoxazole",
        "benzimidazole",
        "indole",
        "benzofuran",
        "quinoline",
        "isoquinoline",
        "benzothiazole",
    ]

    # Track which scaffold is preserved (if any)
    preserved_scaffold = None

    def dfs_traverse(node):
        nonlocal all_mols_have_scaffold, mol_count, preserved_scaffold

        if node["type"] == "mol":
            mol_count += 1
            mol_smiles = node["smiles"]

            # Skip small molecules (likely reagents)
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol and mol.GetNumAtoms() <= 5:
                print(f"Skipping small molecule: {mol_smiles}")
                return

            # Check for heterocyclic scaffolds
            has_scaffold = False

            for scaffold in heterocyclic_scaffolds:
                if checker.check_ring(scaffold, mol_smiles):
                    print(f"Found {scaffold} scaffold in: {mol_smiles}")
                    has_scaffold = True

                    # If this is the first molecule with a scaffold, set it as the preserved one
                    if preserved_scaffold is None:
                        preserved_scaffold = scaffold
                    # If we already found a different scaffold in previous molecules, ensure consistency
                    elif preserved_scaffold != scaffold:
                        print(
                            f"Different scaffold found: expected {preserved_scaffold}, found {scaffold}"
                        )
                        all_mols_have_scaffold = False

                    break

            if not has_scaffold and mol and mol.GetNumAtoms() > 5:
                print(f"Molecule without heterocyclic scaffold: {mol_smiles}")
                all_mols_have_scaffold = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is valid if we have multiple molecules and all have the same scaffold
    strategy_valid = all_mols_have_scaffold and mol_count > 1 and preserved_scaffold is not None

    print(f"Preserved heterocyclic scaffold strategy detected: {strategy_valid}")
    print(f"Number of molecules examined: {mol_count}")
    if preserved_scaffold:
        print(f"Preserved scaffold: {preserved_scaffold}")

    return strategy_valid
