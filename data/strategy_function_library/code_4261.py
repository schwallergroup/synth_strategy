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
    This function detects the use of piperazine as a linker between aromatic fragments.
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

    piperazine_linker_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_linker_detected, findings_json

        if node["type"] == "mol" and node["smiles"]:
            try:
                mol_smiles = node["smiles"]
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    piperazine_indices = checker.get_ring_atom_indices("piperazine", mol_smiles)

                    if piperazine_indices:
                        if "piperazine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("piperazine")

                        for ring_indices in piperazine_indices:
                            piperazine_atoms = set(ring_indices)
                            nitrogen_atoms = []
                            for atom_idx in piperazine_atoms:
                                atom = mol.GetAtomWithIdx(atom_idx)
                                if atom.GetSymbol() == "N":
                                    nitrogen_atoms.append(atom_idx)

                            aromatic_connections = 0
                            for n_idx in nitrogen_atoms:
                                n_atom = mol.GetAtomWithIdx(n_idx)

                                for neighbor in n_atom.GetNeighbors():
                                    neighbor_idx = neighbor.GetIdx()

                                    if neighbor_idx in piperazine_atoms:
                                        continue

                                    if neighbor.GetIsAromatic():
                                        aromatic_connections += 1
                                        if "aromatic_ring" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("aromatic_ring")
                                        break

                            if aromatic_connections >= 2:
                                piperazine_linker_detected = True
                                if {"type": "count", "details": {"target": "aromatic_substituents_on_piperazine", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "aromatic_substituents_on_piperazine", "operator": ">=", "value": 2}})
                                return
            except Exception:
                # Silently ignore errors in SMILES processing
                pass

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is not 'reaction' (e.g., 'chemical' or 'mol')
                new_depth = depth + 1
            # If current node is 'reaction', new_depth remains 'depth'
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return piperazine_linker_detected, findings_json
