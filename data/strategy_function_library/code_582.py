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
    Detects if a brominated pyridine scaffold is preserved throughout the synthesis.
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

    # Track depths with brominated pyridine
    depths_with_bromopyridine = set()
    mol_depths = set()

    def has_bromopyridine(mol_smiles):
        """Check if a molecule contains a brominated pyridine scaffold"""
        nonlocal findings_json
        # First check if it has a pyridine ring
        has_pyridine = checker.check_ring("pyridine", mol_smiles)
        if has_pyridine:
            if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("pyridine")
        else:
            return False

        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Use a detailed approach by finding pyridine ring atoms and checking for attached bromine
        try:
            # Get pyridine ring atoms
            pyridine_indices = checker.get_ring_atom_indices("pyridine", mol_smiles)
            if not pyridine_indices:
                return False

            pyridine_atoms = set()
            # Handle different possible return formats from the checker
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
            found_bromine_on_pyridine = False
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "Br":
                    if "bromine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("bromine")
                    for bond in atom.GetBonds():
                        other_atom_idx = bond.GetOtherAtomIdx(atom.GetIdx())
                        if other_atom_idx in pyridine_atoms:
                            found_bromine_on_pyridine = True
                            break
                if found_bromine_on_pyridine:
                    break

            if found_bromine_on_pyridine:
                # Add the co-occurrence constraint if both are found on the same scaffold
                constraint = {
                    "type": "co-occurrence",
                    "details": {
                        "targets": [
                            "pyridine",
                            "bromine"
                        ],
                        "scope": "same_scaffold"
                    }
                }
                if constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint)
                return True
            return False

        except Exception as e:
            return False

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_bromopyridine, mol_depths, findings_json

        if node["type"] == "mol":
            mol_depths.add(depth)
            if has_bromopyridine(node["smiles"]):
                depths_with_bromopyridine.add(depth)

                # Check if this is a starting material
                if node.get("in_stock", False):
                    pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases only when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if bromopyridine is present at all molecule depths
    is_preserved = len(depths_with_bromopyridine) == len(mol_depths) and len(mol_depths) > 0

    if is_preserved:
        # Add the count constraint if the scaffold is preserved throughout
        constraint = {
            "type": "count",
            "details": {
                "target": "molecules_with_bromopyridine_scaffold",
                "operator": "==",
                "value_target": "total_molecules_in_route"
            }
        }
        if constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint)

    return is_preserved, findings_json
