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
    Detects the reduction of a nitro group to an amine in the early stages of a synthesis, defined as a reaction depth of 3 or greater. The identification is performed using a specific, robust reaction checker for this transformation.
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

    nitro_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected, findings_json
        if nitro_reduction_detected:
            return

        if node["type"] == "reaction" and depth >= 3:  # Early stage (depth >= 3)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                # Direct check using the reaction checker is the most robust method.
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_reduction_detected = True
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Reduction of nitro groups to amines",
                            "position": {
                                "operator": ">=",
                                "value": 3,
                                "unit": "depth"
                            }
                        }
                    })
                    return
            except Exception:
                # Silently ignore reactions that fail to process, e.g., missing rsmi.
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] == 'chemical'
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return nitro_reduction_detected, findings_json
