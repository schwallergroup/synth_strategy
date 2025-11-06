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
    This function detects if the synthesis route maintains a pyrazolopyridine core throughout
    the final product and late-stage intermediates.
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

    all_molecules_have_core = True
    molecule_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal all_molecules_have_core, molecule_count, findings_json

        if node["type"] == "mol" and node.get("smiles"):
            molecule_count += 1

            # Only check for the core in the final product and late-stage intermediates
            # Skip starting materials
            if depth <= 2 and not node.get("in_stock", False):
                # For pyrazolopyridine, we need both rings fused together
                # This SMARTS pattern represents the fused ring system
                pyrazolopyridine_pattern = Chem.MolFromSmarts("c1nn2c(c1)cccc2")
                mol = Chem.MolFromSmiles(node["smiles"])
                has_fused_system = mol.HasSubstructMatch(pyrazolopyridine_pattern) if mol else False

                if has_fused_system:
                    # Record the finding of the pyrazolopyridine ring system
                    if "pyrazolopyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazolopyridine")
                else:
                    all_molecules_have_core = False
                    print(
                        f"Molecule without pyrazolopyridine core at depth {depth}: {node['smiles']}"
                    )

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'mol' or 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Ensure we've actually checked some molecules
    result = all_molecules_have_core and molecule_count > 0

    # If the core was maintained throughout the late stage, add the structural constraint
    if result and "pyrazolopyridine" in findings_json["atomic_checks"]["ring_systems"]:
        # This structural constraint is met if the overall check passes and the ring system was found
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "pyrazolopyridine",
                "position": "late_stage"
            }
        })

    return result, findings_json
