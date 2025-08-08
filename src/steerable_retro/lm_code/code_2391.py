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
    Detects a strategy involving preservation of a pyrazole with difluoromethyl
    substituent throughout the synthesis.
    """
    # Track if we found the motif in the final product and intermediates
    motif_occurrences = []

    # Add depth information to nodes
    def add_depth(node, current_depth=0):
        node["depth"] = current_depth
        for child in node.get("children", []):
            add_depth(child, current_depth + 1)

    # Add depth to all nodes
    add_depth(route)

    def check_pyrazole_with_difluoromethyl(smiles):
        """Helper function to check if a molecule has a pyrazole with difluoromethyl group attached"""
        try:
            # Check for pyrazole
            has_pyrazole = checker.check_ring("pyrazole", smiles)
            if not has_pyrazole:
                return False

            # Create molecule for more detailed analysis
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Get pyrazole ring atoms
            pyrazole_indices = checker.get_ring_atom_indices("pyrazole", smiles)
            if not pyrazole_indices:
                return False

            # Look for a pattern that could be a difluoromethyl group attached to pyrazole
            # This is a more direct approach than relying on functional group detection
            for p_atoms in pyrazole_indices:
                p_atom_set = set(p_atoms)

                # Check each pyrazole atom for a potential difluoromethyl attachment
                for p_atom_idx in p_atoms:
                    p_atom = mol.GetAtomWithIdx(p_atom_idx)

                    # Check neighbors of pyrazole atoms
                    for neighbor in p_atom.GetNeighbors():
                        # If neighbor is carbon, check if it has two fluorine atoms attached
                        if neighbor.GetSymbol() == "C":
                            carbon_idx = neighbor.GetIdx()
                            # Count fluorine neighbors on this carbon
                            f_count = sum(
                                1 for n in neighbor.GetNeighbors() if n.GetSymbol() == "F"
                            )

                            # If carbon has exactly 2 fluorines, it's a difluoromethyl group
                            if f_count == 2:
                                print(f"Found pyrazole with difluoromethyl in {smiles}")
                                return True

            return False
        except Exception as e:
            print(f"Error in check_pyrazole_with_difluoromethyl: {e}")
            return False

    def dfs_traverse(node):
        if node["type"] == "mol":
            smiles = node["smiles"]
            depth = node.get("depth", 0)

            # Check for pyrazole with difluoromethyl
            if check_pyrazole_with_difluoromethyl(smiles):
                print(f"Found pyrazole with difluoromethyl at depth {depth}, SMILES: {smiles}")
                motif_occurrences.append((depth, smiles))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Sort occurrences by depth
    motif_occurrences.sort(key=lambda x: x[0])

    # Check if we have the motif in the final product (depth 0)
    found_in_final = any(depth == 0 for depth, _ in motif_occurrences)

    # Count unique intermediates (depth > 0) with the motif
    intermediate_depths = set(depth for depth, _ in motif_occurrences if depth > 0)
    found_in_intermediates = len(intermediate_depths)

    print(f"Found in final: {found_in_final}, Found in intermediates: {found_in_intermediates}")
    print(f"Total occurrences: {len(motif_occurrences)}")

    # Check if the motif is preserved throughout (found in final and at least 2 intermediates)
    return found_in_final and found_in_intermediates >= 2
