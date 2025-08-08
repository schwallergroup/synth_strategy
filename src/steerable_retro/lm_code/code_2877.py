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
    This function detects a strategy where a trifluoromethyl-substituted isoxazole motif
    is maintained throughout the synthesis.
    """
    # Track presence of the motif at each step
    motif_present_at_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {mol_smiles}")

            try:
                # Check for isoxazole ring
                has_isoxazole = checker.check_ring("isoxazole", mol_smiles)

                # Check for trifluoromethyl group
                has_trifluoro = checker.check_fg("Trifluoro group", mol_smiles)

                if has_isoxazole and has_trifluoro:
                    # Need to verify the trifluoromethyl group is attached to the isoxazole
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if not mol:
                        print(f"Failed to parse SMILES at depth {depth}: {mol_smiles}")
                        return

                    # Get atom indices for isoxazole ring
                    isoxazole_indices_list = checker.get_ring_atom_indices("isoxazole", mol_smiles)

                    # Get atom indices for trifluoromethyl group
                    trifluoro_indices_list = checker.get_fg_atom_indices(
                        "Trifluoro group", mol_smiles
                    )

                    print(f"Isoxazole indices: {isoxazole_indices_list}")
                    print(f"Trifluoromethyl indices: {trifluoro_indices_list}")

                    # Check if any trifluoromethyl group is attached to the isoxazole ring
                    is_attached = False

                    # Flatten isoxazole indices if needed
                    isoxazole_atoms = []
                    if isinstance(isoxazole_indices_list, list):
                        for indices in isoxazole_indices_list:
                            if isinstance(indices, tuple):
                                isoxazole_atoms.extend(indices)
                            else:
                                isoxazole_atoms.append(indices)
                    else:
                        isoxazole_atoms = [isoxazole_indices_list]

                    # Flatten trifluoromethyl indices if needed
                    trifluoro_atoms = []
                    if isinstance(trifluoro_indices_list, list):
                        for indices in trifluoro_indices_list:
                            if isinstance(indices, tuple):
                                trifluoro_atoms.extend(indices)
                            else:
                                trifluoro_atoms.append(indices)
                    else:
                        trifluoro_atoms = [trifluoro_indices_list]

                    # Find the carbon atom of CF3 (usually the first atom in the group)
                    cf3_carbon_atoms = []
                    for i, atom_idx in enumerate(trifluoro_atoms):
                        if isinstance(atom_idx, int):
                            atom = mol.GetAtomWithIdx(atom_idx)
                            if atom.GetSymbol() == "C":
                                cf3_carbon_atoms.append(atom_idx)

                    # Check if any CF3 carbon is connected to any isoxazole atom
                    for cf3_carbon in cf3_carbon_atoms:
                        carbon_atom = mol.GetAtomWithIdx(cf3_carbon)
                        for neighbor in carbon_atom.GetNeighbors():
                            neighbor_idx = neighbor.GetIdx()
                            if neighbor_idx in isoxazole_atoms:
                                is_attached = True
                                break

                        # Also check if any isoxazole atom is connected to the CF3 carbon
                        if not is_attached:
                            for isox_idx in isoxazole_atoms:
                                isox_atom = mol.GetAtomWithIdx(isox_idx)
                                for neighbor in isox_atom.GetNeighbors():
                                    if neighbor.GetIdx() == cf3_carbon:
                                        is_attached = True
                                        break
                                if is_attached:
                                    break

                        if is_attached:
                            break

                    if is_attached:
                        motif_present_at_steps.append(depth)
                        print(f"Detected trifluoromethyl-isoxazole at depth {depth}")
                    else:
                        print(
                            f"Found isoxazole and trifluoromethyl, but they are not connected at depth {depth}"
                        )
                else:
                    if not has_isoxazole:
                        print(f"No isoxazole ring found at depth {depth}")
                    if not has_trifluoro:
                        print(f"No trifluoromethyl group found at depth {depth}")
            except Exception as e:
                print(f"Error processing molecule at depth {depth}: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    print(f"Motif found at depths: {motif_present_at_steps}")

    # Strategy is present if the motif appears in at least 3 different depths
    # and spans from early to late stages of the synthesis
    if len(motif_present_at_steps) >= 3:
        depth_range = max(motif_present_at_steps) - min(motif_present_at_steps)
        print(f"Depth range: {depth_range}")
        if depth_range >= 2:  # Spans at least 3 steps
            print("Strategy detected: Trifluoromethyl-isoxazole maintained throughout synthesis")
            return True
        else:
            print("Motif doesn't span enough synthesis steps")
    else:
        print(
            f"Motif not found in enough steps (found in {len(motif_present_at_steps)} steps, need at least 3)"
        )

    return False
