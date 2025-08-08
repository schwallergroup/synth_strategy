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
    Detects a synthetic strategy involving nitro reduction.
    Looks for conversion of aromatic nitro group to amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {depth}")
                    has_nitro_reduction = True
                else:
                    # Alternative check: look for nitro groups in reactants and amines in product
                    reactant_with_nitro = None
                    for r in reactants_smiles:
                        if checker.check_fg("Nitro group", r):
                            reactant_with_nitro = r
                            break

                    if reactant_with_nitro and checker.check_fg("Primary amine", product_smiles):
                        # Check if the nitro group is on an aromatic ring
                        reactant_mol = Chem.MolFromSmiles(reactant_with_nitro)
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if reactant_mol and product_mol:
                            # Method 1: Check using atom indices
                            try:
                                nitro_indices = (
                                    checker.get_fg_atom_indices("Nitro group", reactant_with_nitro)
                                    or []
                                )

                                for atom_indices_group in nitro_indices:
                                    for atom_indices in atom_indices_group:
                                        # Check if we have valid atom indices
                                        if (
                                            atom_indices
                                            and isinstance(atom_indices, tuple)
                                            and len(atom_indices) > 0
                                        ):
                                            # Get the atom attached to the nitro group
                                            try:
                                                n_atom = reactant_mol.GetAtomWithIdx(
                                                    atom_indices[0]
                                                )
                                                for neighbor in n_atom.GetNeighbors():
                                                    if neighbor.GetIsAromatic():
                                                        print(
                                                            f"Found aromatic nitro reduction at depth {depth}"
                                                        )
                                                        has_nitro_reduction = True
                                                        break
                                            except Exception as e:
                                                print(f"Error processing atom indices: {e}")
                            except Exception as e:
                                print(f"Error with atom indices approach: {e}")

                            # Method 2: Direct check for aromatic nitro groups
                            if not has_nitro_reduction:
                                for atom in reactant_mol.GetAtoms():
                                    if atom.GetIsAromatic():
                                        for neighbor in atom.GetNeighbors():
                                            if (
                                                neighbor.GetSymbol() == "N"
                                                and neighbor.GetIdx()
                                                in [
                                                    a.GetIdx()
                                                    for a in reactant_mol.GetAtoms()
                                                    if a.GetSymbol() == "N"
                                                ]
                                            ):
                                                # Check if this nitrogen is part of a nitro group
                                                oxygen_count = sum(
                                                    1
                                                    for nn in neighbor.GetNeighbors()
                                                    if nn.GetSymbol() == "O"
                                                )
                                                if (
                                                    oxygen_count >= 2
                                                ):  # Nitro group typically has 2 oxygen atoms
                                                    print(
                                                        f"Found aromatic nitro reduction at depth {depth} (direct check)"
                                                    )
                                                    has_nitro_reduction = True
                                                    break
                                        if has_nitro_reduction:
                                            break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if nitro reduction is found
    return has_nitro_reduction
