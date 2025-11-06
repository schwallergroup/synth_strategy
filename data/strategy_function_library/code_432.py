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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid to amide conversion",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route follows a linear strategy with
    a late-stage amide bond disconnection.
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

    late_stage_disconnection = False
    is_linear = True
    total_nodes = 0
    branching_nodes = 0

    # First pass to check if the route is mostly linear
    def check_linearity(node):
        nonlocal total_nodes, branching_nodes

        total_nodes += 1
        if len(node.get("children", [])) > 1:
            branching_nodes += 1

        for child in node.get("children", []):
            check_linearity(child)

    # Second pass to find late-stage amide disconnection
    def dfs_traverse(node, depth=0):
        nonlocal late_stage_disconnection, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check for amide bond disconnection at low depth (late stage)
            if depth <= 3:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                is_amide_coupling = False
                for name in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_amide_coupling:
                    late_stage_disconnection = True
                    # Add positional constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "amide_coupling",
                            "position": "depth <= 3"
                        }
                    })

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        # Process children nodes
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Check if the route is mostly linear
    check_linearity(route)
    linearity_ratio = 1.0 - (branching_nodes / max(1, total_nodes))
    
    # Add linearity constraint if met
    if linearity_ratio >= 0.7:
        is_linear = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "route_linearity_ratio",
                "operator": ">=",
                "value": 0.7
            }
        })
    else:
        is_linear = False

    # Only proceed if the route is mostly linear
    if is_linear:
        dfs_traverse(route)

    result = is_linear and late_stage_disconnection
    
    # Add co-occurrence constraint if met
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "route_linearity_ratio >= 0.7",
                    "late_stage_amide_coupling"
                ]
            }
        })

    return result, findings_json
