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


N_METHYLATION_REACTIONS = [
    "N-methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "Methylation with MeI_primary",
    "Methylation with MeI_secondary",
    "Methylation with MeI_tertiary",
    "DMS Amine methylation",
    "Parnes methylation",
    "Methylation with DMS",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where specific N-methylation reactions are performed in the late stages of a synthesis. The strategy is flagged if at least two such reactions occur, with at least 50% of them happening in the latter half of the synthesis (i.e., at low retrosynthetic depth). The function identifies N-methylations by checking for a predefined list of reaction types, including Eschweiler-Clarke, reductive amination with formaldehyde, and methylation with common reagents like MeI and DMS.
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

    # Track N-methylations at different depths
    n_methylations_by_depth = {}
    max_depth = 0

    def is_n_methylation(rsmi):
        """Checks if a reaction corresponds to a known N-methylation reaction type from a predefined list."""
        for reaction_type in N_METHYLATION_REACTIONS:
            if checker.check_reaction(reaction_type, rsmi):
                if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal n_methylations_by_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                if is_n_methylation(rsmi):
                    if depth not in n_methylations_by_depth:
                        n_methylations_by_depth[depth] = 0
                    n_methylations_by_depth[depth] += 1
            except KeyError:
                pass
            except Exception:
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Calculate the percentage of N-methylations in the first half of the synthesis
    total_methylations = sum(n_methylations_by_depth.values())
    if total_methylations == 0:
        return False, findings_json

    # Define "late stage" as the first half of the synthesis (lower depths)
    late_stage_threshold = max_depth / 2
    late_stage_methylations = sum(
        count for depth, count in n_methylations_by_depth.items() if depth <= late_stage_threshold
    )

    # Strategy criteria: at least 2 N-methylations with at least 50% in late stage
    late_stage_percentage = (
        late_stage_methylations / total_methylations if total_methylations > 0 else 0
    )
    strategy_detected = total_methylations >= 2 and late_stage_percentage >= 0.5

    if total_methylations >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "N-methylation",
                "operator": ">=",
                "value": 2
            }
        })
    
    if late_stage_percentage >= 0.5:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "N-methylation",
                "qualifier": "ratio",
                "position_definition": "depth <= max_depth / 2",
                "operator": ">=",
                "value": 0.5
            }
        })

    return strategy_detected, findings_json
