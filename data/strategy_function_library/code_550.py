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
    This function detects a late-stage nitro reduction strategy where a nitro group
    is reduced to an amine in the final or near-final step of the synthesis.
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

    has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction, findings_json

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            try:
                # Check if this is a nitro reduction reaction using the specific checker
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                if is_nitro_reduction:
                    print(f"Found nitro reduction reaction at depth {depth}")
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                    # Check if this is a late-stage reaction (depth 1 or 2)
                    if depth <= 2:
                        print(f"Found late-stage nitro reduction at depth {depth}")
                        has_nitro_reduction = True
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Reduction of nitro groups to amines",
                                "position": "late_stage",
                                "description": "The reaction must occur at a retrosynthesis depth of 2 or less."
                            }
                        })
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child_idx, child in enumerate(node.get("children", [])):
            # Increment depth based on the new rule
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node
                next_depth = depth
            else:
                # Depth increases when traversing from a chemical node
                next_depth = depth + 1
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {has_nitro_reduction}")

    return has_nitro_reduction, findings_json
