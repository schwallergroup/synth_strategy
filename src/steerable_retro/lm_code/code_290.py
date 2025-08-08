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
    This function detects a synthetic strategy involving imine formation with a highly branched amine.
    """
    # Initialize tracking variables
    has_branched_amine_imine = False

    def is_highly_branched_amine(mol_smiles):
        """Check if the molecule contains a highly branched primary amine"""
        if not checker.check_fg("Primary amine", mol_smiles):
            return False

        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Get primary amine nitrogen atoms
        amine_indices = []
        try:
            for atom_indices_group in checker.get_fg_atom_indices("Primary amine", mol_smiles):
                if isinstance(atom_indices_group, tuple) and all(
                    isinstance(x, tuple) for x in atom_indices_group
                ):
                    # Handle case where we get a tuple of tuples
                    for indices in atom_indices_group:
                        if indices and len(indices) > 0:
                            amine_indices.append(indices[0])
                elif isinstance(atom_indices_group, tuple) and len(atom_indices_group) > 0:
                    # Handle case where we get a single tuple of indices
                    amine_indices.append(atom_indices_group[0])
                elif isinstance(atom_indices_group, int):
                    # Handle case where we get a single integer
                    amine_indices.append(atom_indices_group)
        except Exception as e:
            print(f"Error processing amine indices: {e}")
            # Try a simpler approach - just find nitrogen atoms with one hydrogen
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 7 and atom.GetTotalNumHs() == 1:
                    amine_indices.append(atom.GetIdx())

        if not amine_indices:
            print(f"No primary amine nitrogens found in {mol_smiles}")
            return False

        print(f"Found {len(amine_indices)} primary amine nitrogens in {mol_smiles}")

        for n_idx in amine_indices:
            n_atom = mol.GetAtomWithIdx(n_idx)
            # Get the carbon attached to the nitrogen
            for neighbor in n_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    carbon_idx = neighbor.GetIdx()
                    carbon = mol.GetAtomWithIdx(carbon_idx)

                    # Count non-hydrogen branches from this carbon
                    branch_count = 0
                    for c_neighbor in carbon.GetNeighbors():
                        if (
                            c_neighbor.GetAtomicNum() != 1 and c_neighbor.GetIdx() != n_idx
                        ):  # Not H and not the N
                            branch_count += 1

                            # Check if any of these branches are themselves branched
                            sub_branch_count = 0
                            for branch_neighbor in c_neighbor.GetNeighbors():
                                if (
                                    branch_neighbor.GetAtomicNum() != 1
                                    and branch_neighbor.GetIdx() != carbon_idx
                                ):
                                    sub_branch_count += 1

                            if sub_branch_count > 1:  # If the branch has multiple non-H attachments
                                branch_count += 1  # Count it as an additional branch

                    if branch_count >= 2:  # Consider it highly branched if it has 2+ branches
                        print(
                            f"Found highly branched amine with {branch_count} branches in {mol_smiles}"
                        )
                        return True

        return False

    def dfs_traverse(node):
        nonlocal has_branched_amine_imine

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check for imine formation reaction
                imine_formation = checker.check_reaction(
                    "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                ) or checker.check_reaction(
                    "Addition of primary amines to ketones/thiocarbonyls", rsmi
                )

                if imine_formation:
                    print(f"Detected imine formation reaction: {rsmi}")

                # Verify product has imine and reactants don't
                has_imine_product = checker.check_fg(
                    "Substituted imine", product_smiles
                ) or checker.check_fg("Unsubstituted imine", product_smiles)
                no_imine_reactants = not any(
                    checker.check_fg("Substituted imine", r)
                    or checker.check_fg("Unsubstituted imine", r)
                    for r in reactants_smiles
                )

                if has_imine_product:
                    print(f"Product has imine: {product_smiles}")

                if imine_formation or (has_imine_product and no_imine_reactants):
                    # Check for branched amine reactant
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amine", r):
                            print(f"Checking if reactant is a highly branched amine: {r}")
                            if is_highly_branched_amine(r):
                                print(f"Detected imine formation with branched amine")
                                print(f"Reaction: {rsmi}")
                                has_branched_amine_imine = True
                                break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Strategy detection results: has_branched_amine_imine={has_branched_amine_imine}")

    return has_branched_amine_imine
