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
    Detects a late-stage N-acylation, specifically the reaction of a nitrogen nucleophile
    with an acyl halide, thioacyl halide, carbamoyl halide, or a related acylating agent.
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

    has_late_stage_n_acylation = False

    def is_n_acylation_with_acyl_halide(reaction_node, depth):
        """Check if a reaction is an N-acylation with acyl halide"""
        nonlocal findings_json
        if "rsmi" not in reaction_node.get("metadata", {}):
            return False

        # Only consider reactions at depth 0 or 1 (late stage)
        if depth > 1:
            return False

        rsmi = reaction_node["metadata"]["mapped_reaction_smiles"]

        # Check if this is an acylation reaction using the checker function
        target_reaction_name = "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N"
        if checker.check_reaction(
            target_reaction_name,
            rsmi,
        ):
            print(f"Found N-acylation reaction at depth {depth}")
            findings_json["atomic_checks"]["named_reactions"].append(target_reaction_name)
            # Add structural constraint if it's late stage
            if depth <= 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": target_reaction_name,
                        "position": "late_stage"
                    }
                })
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_n_acylation, findings_json

        if node["type"] == "reaction":
            if is_n_acylation_with_acyl_halide(node, depth):
                has_late_stage_n_acylation = True
                print(
                    f"Found late-stage N-acylation at depth {depth}: {node.get('metadata', {}).get('rsmi', 'No RSMI')}"
                )

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, increase depth
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Late-stage N-acylation found: {has_late_stage_n_acylation}")

    # We keep the original function name for compatibility
    return has_late_stage_n_acylation, findings_json
