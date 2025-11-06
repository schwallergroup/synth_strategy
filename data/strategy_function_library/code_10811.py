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


LATE_STAGE_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Sonogashira alkyne_aryl halide",
    "Buchwald-Hartwig",
    "Heck_terminal_vinyl",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage, multi-component coupling reaction, where three or more significant fragments (>= 5 heavy atoms) are joined in the final synthetic step (depth=1). The reaction must be one of the specified named coupling types defined in LATE_STAGE_COUPLING_REACTIONS.
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

    found_late_coupling = False

    def count_significant_fragments(reactants_smiles):
        """Count reactants that are significant fragments (not small reagents)"""
        significant_fragments = 0
        for smiles in reactants_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.GetNumHeavyAtoms() >= 5:  # Consider as significant if 5+ heavy atoms
                significant_fragments += 1
        return significant_fragments

    def infer_depth(node, current_depth=0):
        """Infer depth if not provided in metadata"""
        if node["type"] == "mol" and node.get("in_stock", False):
            return current_depth

        max_child_depth = current_depth
        for child in node.get("children", []):
            child_depth = infer_depth(child, current_depth + 1)
            max_child_depth = max(max_child_depth, child_depth)

        return max_child_depth

    def is_coupling_reaction(rsmi):
        """Check if the reaction is one of the predefined named coupling reactions."""
        try:
            for reaction_name in LATE_STAGE_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    return True
            return False
        except:
            return False

    def dfs_traverse(node, calculated_depth=None):
        nonlocal found_late_coupling, findings_json

        # Determine the depth for the current node based on its type and the calculated_depth
        current_node_depth = calculated_depth
        if node["type"] == "reaction":
            # If it's a reaction node, use its metadata depth if available, otherwise the calculated depth
            metadata_depth = node.get("metadata", {}).get("depth", -1)
            if metadata_depth != -1:
                current_node_depth = metadata_depth

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Check for late-stage fragment coupling (depth 1 is the final step)
            if current_node_depth == 1:
                # Record positional constraint if depth is 1
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "named_coupling_reaction",
                        "position": "last_stage"
                    }
                })

                # Record count constraint for reactants
                if len(reactants) >= 3:
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants",
                            "operator": ">=",
                            "value": 3
                        }
                    })

                # Count significant fragments (not small reagents)
                significant_count = count_significant_fragments(reactants)

                if significant_count >= 3:
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants_with_gte_5_heavy_atoms",
                            "operator": ">=",
                            "value": 3
                        }
                    })

                if significant_count >= 3 and is_coupling_reaction(rsmi):
                    found_late_coupling = True
                    print(
                        f"Found late-stage coupling with {significant_count} significant fragments at depth {current_node_depth}"
                    )

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            next_depth = current_node_depth
            if node["type"] != "reaction": # This means current node is 'chemical' (or 'mol')
                next_depth += 1
            
            dfs_traverse(child, next_depth)

    # Calculate max depth for the route if needed
    max_depth = infer_depth(route)

    # Start traversal with calculated depth
    dfs_traverse(route, 0)

    print(f"Late-stage fragment coupling strategy detected: {found_late_coupling}")
    return found_late_coupling, findings_json
