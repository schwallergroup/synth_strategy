from typing import Tuple, Dict, List
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
    This function detects a linear synthesis strategy that maintains a trifluoromethyl-pyrazole
    moiety throughout the synthesis.
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

    # Track if we found the CF3-pyrazole in all steps
    steps_with_cf3_pyrazole = set()
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_steps, findings_json

        node_type = node.get("type")

        if node_type == "reaction":
            # Count all reaction nodes
            total_steps += 1

            # Get reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Extract product from reaction SMILES
            try:
                parts = rsmi.split(">")
                if len(parts) < 3:
                    return
                product = parts[2]
            except Exception as e:
                return

            # Get depth for tracking
            # The depth passed to the recursive call should be `depth` if the current node's type is 'reaction',
            # and `depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
            # The current node's depth is `depth` as passed to this function.
            node["metadata"] = node.get("metadata", {})
            node["metadata"]["depth"] = depth

            # Check if product contains both trifluoromethyl group and pyrazole ring
            has_cf3 = checker.check_fg("Trifluoro group", product)
            has_pyrazole = checker.check_ring("pyrazole", product)

            if has_cf3:
                if "Trifluoro group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Trifluoro group")
            if has_pyrazole:
                if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

            if has_cf3 and has_pyrazole:
                # Verify they are connected in the molecule
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Get atom indices for both substructures
                    cf3_indices = checker.get_fg_atom_indices("Trifluoro group", product)
                    pyrazole_indices = checker.get_ring_atom_indices("pyrazole", product)

                    # Check if any CF3 group is connected to a pyrazole ring
                    if cf3_indices and pyrazole_indices:
                        connected = False

                        # Iterate through all CF3 groups
                        for cf3_group in cf3_indices:
                            # Extract the atom indices for this CF3 group
                            cf3_atoms = set(cf3_group)

                            # Iterate through all pyrazole rings
                            for pyrazole_ring in pyrazole_indices:
                                # Extract the atom indices for this pyrazole ring
                                pyrazole_atoms = set(pyrazole_ring)

                                # Check for bonds between atoms (adjacent connection)
                                for c in cf3_atoms:
                                    for p in pyrazole_atoms:
                                        bond = product_mol.GetBondBetweenAtoms(c, p)
                                        if bond is not None:
                                            connected = True
                                            break
                                    if connected:
                                        break

                            if connected:
                                break

                        if connected:
                            steps_with_cf3_pyrazole.add(depth)

        # Traverse children
        for child in node.get("children", []):
            new_depth = depth
            if node_type != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if CF3-pyrazole is present in all steps
    result = len(steps_with_cf3_pyrazole) == total_steps and total_steps > 0

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_product_lacking_connected_trifluoromethyl_pyrazole",
                "operator": "==",
                "value": 0
            }
        })

    return result, findings_json
