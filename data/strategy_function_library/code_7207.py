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
    """Detects syntheses that involve the formation of pyrazole or imidazole rings. The strategy is flagged if these rings are formed, with a special focus on formations occurring in the early stages of the synthesis."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Track heterocycle occurrences at different depths
    pyrazole_depths = []
    imidazole_depths = []
    total_nodes = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_nodes, pyrazole_depths, imidazole_depths, findings_json
        total_nodes += 1

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Check if the reaction involves formation of heterocycles
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if pyrazole is formed in the reaction
            if not checker.check_ring("pyrazole", reactants) and checker.check_ring(
                "pyrazole", product
            ):
                pyrazole_depths.append(depth)
                if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Check if imidazole is formed in the reaction
            if not checker.check_ring("imidazole", reactants) and checker.check_ring(
                "imidazole", product
            ):
                imidazole_depths.append(depth)
                if "imidazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("imidazole")
                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Analyze results
    pyrazole_count = len(pyrazole_depths)
    imidazole_count = len(imidazole_depths)

    print(f"Total nodes: {total_nodes}")
    print(f"Pyrazole occurrences: {pyrazole_count} at depths {pyrazole_depths}")
    print(f"Imidazole occurrences: {imidazole_count} at depths {imidazole_depths}")

    # Determine if the synthesis is heterocycle-based
    # Criteria:
    # 1. At least one heterocycle is present
    # 2. Heterocycles appear in early stages (higher depths) or are formed in reactions
    # 3. Heterocycles represent a significant portion of the synthesis

    result = False

    if pyrazole_count + imidazole_count == 0:
        return result, findings_json

    # Check if heterocycles appear in early stages (higher depths)
    early_stage_threshold = max(
        2, total_nodes // 3
    )  # Consider depths > 2 or top third as early stage
    early_stage_heterocycles = sum(
        1 for d in pyrazole_depths + imidazole_depths if d >= early_stage_threshold
    )

    # Check if heterocycles represent a significant portion
    heterocycle_ratio = (pyrazole_count + imidazole_count) / max(1, total_nodes)

    print(f"Early stage heterocycles: {early_stage_heterocycles}")
    print(f"Heterocycle ratio: {heterocycle_ratio:.2f}")

    # Return True if either:
    # - Many heterocycles throughout the synthesis (>20% of nodes)
    # - Some heterocycles in early stages (building blocks)
    if heterocycle_ratio > 0.2:
        result = True
        findings_json["structural_constraints"].append({
            "type": "disjunction",
            "details": {
                "conditions": [
                    {
                        "type": "count",
                        "details": {
                            "target": "heterocycle_formation_ratio",
                            "operator": ">",
                            "value": 0.2
                        }
                    }
                ]
            }
        })
    
    if early_stage_heterocycles >= 1:
        result = True
        # If the first condition already added the disjunction, we need to ensure we don't duplicate it
        # and instead add the second condition to the existing disjunction or create a new one if not present.
        # For simplicity, we'll just add the full disjunction if it's not already there, assuming the problem implies
        # that if either condition is met, the *entire* disjunction constraint is considered met and recorded.
        # A more robust solution would check for the existence of the disjunction and append to its 'conditions' list.
        # Given the prompt's example, we'll append the full constraint if this condition is met.
        # To avoid duplicate entries for the *same* constraint object, we'll check if it's already present.
        positional_constraint = {
            "type": "positional",
            "details": {
                "target": "heterocycle_formation",
                "position": "early_stage",
                "min_count": 1
            }
        }
        full_disjunction_constraint = {
            "type": "disjunction",
            "details": {
                "conditions": [
                    {
                        "type": "count",
                        "details": {
                            "target": "heterocycle_formation_ratio",
                            "operator": ">",
                            "value": 0.2
                        }
                    },
                    positional_constraint
                ]
            }
        }

        # Check if the disjunction (or its components) is already recorded to avoid duplicates
        # This is a simplified check; a real-world scenario might require more sophisticated comparison
        # of dictionary contents.
        if not any(d.get('type') == 'disjunction' and 
                   any(c.get('type') == 'positional' and c.get('details', {}).get('position') == 'early_stage' 
                       for c in d.get('details', {}).get('conditions', []))
                   for d in findings_json["structural_constraints"]):
            # If the positional part of the disjunction is not yet recorded, add the full disjunction.
            # This assumes that if one part of the disjunction is met, the whole disjunction is recorded.
            # If the ratio condition was met first, the disjunction would have been added. We need to ensure
            # the positional part is also reflected if it's met.
            # A more precise implementation would be to add only the specific met condition to the 'conditions' list
            # of the disjunction if the disjunction already exists.
            # For this problem, we'll add the full disjunction if the positional condition is met and the disjunction
            # (with the positional part) isn't already there.
            if not any(c == full_disjunction_constraint for c in findings_json["structural_constraints"]):
                findings_json["structural_constraints"].append(full_disjunction_constraint)

    return result, findings_json