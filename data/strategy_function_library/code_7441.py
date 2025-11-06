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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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
    This function detects if the synthesis involves chloropyridazine
    intermediates as key building blocks.
    """
    from rdkit import Chem
    chloropyridazine_present = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth, is_final_product=True):
        nonlocal chloropyridazine_present, findings_json

        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check if this is an intermediate (not a starting material and not the final product)
            is_intermediate = not node.get("in_stock", False) and not is_final_product

            if is_intermediate:
                # Check if the molecule contains a pyridazine ring
                has_pyridazine = checker.check_ring("pyridazine", smiles)

                if has_pyridazine:
                    if "pyridazine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridazine")

                    # Check if it has an aromatic halide
                    has_aromatic_halide = checker.check_fg("Aromatic halide", smiles)

                    if has_aromatic_halide:
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                        # Verify it's specifically a chlorine on the pyridazine
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            # Find pyridazine atoms
                            pyridazine_indices = checker.get_ring_atom_indices("pyridazine", smiles)

                            # Find aromatic halide atoms
                            aromatic_halide_indices = checker.get_fg_atom_indices(
                                "Aromatic halide", smiles
                            )

                            # Check if any aromatic halide is attached to the pyridazine ring
                            for halide_group in aromatic_halide_indices:
                                aromatic_atom = halide_group[0]  # The aromatic carbon
                                halide_atom = halide_group[1]  # The halide atom

                                # Check if the aromatic atom is part of pyridazine
                                for pyridazine_group in pyridazine_indices:
                                    # Check if the aromatic atom is in the pyridazine ring
                                    if aromatic_atom in pyridazine_group:
                                        # Check if the halide is chlorine
                                        atom = mol.GetAtomWithIdx(halide_atom)
                                        if atom.GetSymbol() == "Cl":
                                            chloropyridazine_present = True
                                            # Add the structural constraint if detected
                                            if {
                                                "type": "co-occurrence",
                                                "details": {
                                                    "targets": [
                                                        "pyridazine",
                                                        "Aromatic halide"
                                                    ],
                                                    "scope": "single_intermediate_molecule",
                                                    "condition": "The 'Aromatic halide' must be a chlorine atom attached to the 'pyridazine' ring."
                                                }
                                            } not in findings_json["structural_constraints"]:
                                                findings_json["structural_constraints"].append({
                                                    "type": "co-occurrence",
                                                    "details": {
                                                        "targets": [
                                                            "pyridazine",
                                                            "Aromatic halide"
                                                        ],
                                                        "scope": "single_intermediate_molecule",
                                                        "condition": "The 'Aromatic halide' must be a chlorine atom attached to the 'pyridazine' ring."
                                                    }
                                                })
                                            break # Found a chloropyridazine, no need to check other halide groups
                                    if chloropyridazine_present: # If found, break from inner loop too
                                        break
                                if chloropyridazine_present: # If found, break from outer loop too
                                    break

        # Traverse children (all children are not final products)
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth, is_final_product=False)

    # Start traversal
    dfs_traverse(route, 0) # Initial depth is 0

    return chloropyridazine_present, findings_json
