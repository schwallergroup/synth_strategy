from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the strategy of synthesizing a nitrated pyrazole derivative.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    has_pyrazole = False
    has_nitration = False
    nitration_depth = -1
    nitrated_pyrazole = False
    max_depth_seen = 0  # Track maximum depth to understand early vs late stage

    def is_nitro_attached_to_pyrazole(mol_smiles):
        """Helper function to check if a nitro group is attached to a pyrazole ring"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                print(f"Could not parse SMILES: {mol_smiles}")
                return False

            # Get atom indices for pyrazole rings and nitro groups
            pyrazole_indices = checker.get_ring_atom_indices("pyrazole", mol_smiles)
            nitro_indices = checker.get_fg_atom_indices("Nitro group", mol_smiles)

            if not pyrazole_indices or not nitro_indices:
                return False

            # Convert indices to sets for easier checking
            pyrazole_atoms = set()
            for ring_indices in pyrazole_indices:
                # Handle different return formats
                if isinstance(ring_indices, tuple) and len(ring_indices) > 0:
                    if isinstance(ring_indices[0], tuple):
                        for atom_idx in ring_indices[0]:
                            pyrazole_atoms.add(atom_idx)
                    else:
                        for atom_idx in ring_indices:
                            pyrazole_atoms.add(atom_idx)
                elif isinstance(ring_indices, int):
                    pyrazole_atoms.add(ring_indices)

            # Check if any nitro group is attached to the pyrazole ring
            for nitro_group in nitro_indices:
                nitro_atoms = set()
                # Handle different return formats
                if isinstance(nitro_group, tuple) and len(nitro_group) > 0:
                    if isinstance(nitro_group[0], tuple):
                        for atom_idx in nitro_group[0]:
                            nitro_atoms.add(atom_idx)
                    else:
                        for atom_idx in nitro_group:
                            nitro_atoms.add(atom_idx)
                elif isinstance(nitro_group, int):
                    nitro_atoms.add(nitro_group)

                # Check for adjacency between nitro and pyrazole atoms
                for nitro_atom_idx in nitro_atoms:
                    atom = mol.GetAtomWithIdx(nitro_atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in pyrazole_atoms:
                            return True

                # Also check if pyrazole atoms have nitro neighbors
                for pyrazole_atom_idx in pyrazole_atoms:
                    atom = mol.GetAtomWithIdx(pyrazole_atom_idx)
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() in nitro_atoms:
                            return True

            return False
        except Exception as e:
            print(f"Error checking nitro-pyrazole attachment: {e}")
            return False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_pyrazole, has_nitration, nitration_depth, nitrated_pyrazole, max_depth_seen, findings_json

        # Update max depth seen
        max_depth_seen = max(max_depth_seen, current_depth)

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                # Continue traversal even if rsmi is missing
                for child in node.get("children", []):
                    dfs_traverse(child, current_depth + 1) # Depth increases from reaction to chemical
                return

            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyrazole in product
                product_has_pyrazole = checker.check_ring("pyrazole", product)
                if product_has_pyrazole:
                    has_pyrazole = True
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                    print(f"Found pyrazole at depth {current_depth}")

                # Check for nitro group in product
                product_has_nitro = checker.check_fg("Nitro group", product)
                if product_has_nitro:
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")

                # Check if product has nitrated pyrazole
                if (
                    product_has_pyrazole
                    and product_has_nitro
                    and is_nitro_attached_to_pyrazole(product)
                ):
                    nitrated_pyrazole = True
                    print(f"Found nitrated pyrazole at depth {current_depth}")

                # Check for explicit nitration reaction
                explicit_nitration_reactions = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Non-aromatic nitration with HNO3"
                ]
                is_explicit_nitration = False
                for r_name in explicit_nitration_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        is_explicit_nitration = True
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                # Record nitration if found
                if is_explicit_nitration:
                    has_nitration = True
                    nitration_depth = current_depth
                    print(f"Found nitration reaction at depth {current_depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction
            # If current node is 'reaction', depth remains the same for its children (which are chemicals)
            # If current node is 'chemical', depth increases for its children (which are reactions)
            next_depth = current_depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                next_depth = current_depth + 1
            
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if late-stage nitration (relative to max depth)
    is_late_stage = False
    if nitration_depth >= 0:
        # If max_depth is small, any nitration is late-stage
        if max_depth_seen <= 2:
            is_late_stage = True
        # Otherwise, nitration should be in first third of synthesis
        else:
            is_late_stage = nitration_depth <= (max_depth_seen // 3)

    # Check if we found a nitrated pyrazole with nitration at late stage
    strategy_detected = (
        has_pyrazole
        and nitrated_pyrazole
        and (
            (has_nitration and is_late_stage)  # Explicit late-stage nitration
            or (
                not has_nitration and nitrated_pyrazole
            )  # Nitrated pyrazole without explicit nitration
        )
    )

    if strategy_detected:
        print(
            f"Detected nitrated pyrazole strategy with {'late-stage' if has_nitration else 'implicit'} nitration"
        )
        # Add structural constraints if strategy is detected
        if nitrated_pyrazole:
            # Co-occurrence constraint
            if {
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "pyrazole",
                        "Nitro group"
                    ],
                    "comment": "The function specifically checks that a nitro group is attached to a pyrazole ring in a product molecule."
                }
            } not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "pyrazole",
                            "Nitro group"
                        ],
                        "comment": "The function specifically checks that a nitro group is attached to a pyrazole ring in a product molecule."
                    }
                })
        
        if has_nitration and is_late_stage:
            # Positional constraint
            if {
                "type": "positional",
                "details": {
                    "target": [
                        "Aromatic nitration with HNO3",
                        "Aromatic nitration with NO3 salt",
                        "Aromatic nitration with NO2 salt",
                        "Aromatic nitration with alkyl NO2",
                        "Non-aromatic nitration with HNO3"
                    ],
                    "position": "late_stage",
                    "comment": "If any explicit nitration reaction is found, it must occur late in the synthesis (defined as within the first third of reaction steps from the final product)."
                }
            } not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": [
                            "Aromatic nitration with HNO3",
                            "Aromatic nitration with NO3 salt",
                            "Aromatic nitration with NO2 salt",
                            "Aromatic nitration with alkyl NO2",
                            "Non-aromatic nitration with HNO3"
                        ],
                        "position": "late_stage",
                        "comment": "If any explicit nitration reaction is found, it must occur late in the synthesis (defined as within the first third of reaction steps from the final product)."
                    }
                })

    else:
        print(
            f"Strategy not detected: has_pyrazole={has_pyrazole}, has_nitration={has_nitration}, nitrated_pyrazole={nitrated_pyrazole}, nitration_depth={nitration_depth}, max_depth={max_depth_seen}"
        )

    return strategy_detected, findings_json
