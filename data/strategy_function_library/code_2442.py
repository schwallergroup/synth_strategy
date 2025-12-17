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


LATE_STAGE_AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies late-stage (depth <= 1) amide coupling reactions by matching against the specific named reactions defined in the LATE_STAGE_AMIDE_COUPLING_REACTIONS list.
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

    found_late_amide = False

    # Define the structural constraint object from the strategy JSON
    # This is hardcoded based on the provided strategy JSON for this specific problem.
    positional_constraint = {
      "type": "positional",
      "details": {
        "targets": [
          "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
          "Carboxylic acid with primary amine to amide",
          "Ester with primary amine to amide",
          "Ester with secondary amine to amide",
          "Acyl chloride with primary amine to amide (Schotten-Baumann)",
          "Acyl chloride with secondary amine to amide",
          "Acylation of primary amines",
          "Acylation of secondary amines"
        ],
        "logic": "ANY",
        "position": {
          "variable": "depth",
          "operator": "<=",
          "value": 1
        }
      }
    }

    def dfs_traverse(node, current_depth=0):
        nonlocal found_late_amide, findings_json

        # Store depth for future reference
        if "depth" not in node:
            node["depth"] = current_depth

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            depth = node.get("depth", current_depth)

            if depth <= 1:  # Late stage (depth 0 or 1)
                # Check for amide coupling reactions using the checker
                for reaction_type in LATE_STAGE_AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        found_late_amide = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add the structural constraint if a late-stage amide coupling is found
                        if positional_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(positional_constraint)
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                next_depth = current_depth  # Depth remains the same from reaction to chemical
            else:
                next_depth = current_depth + 1 # Depth increases from chemical to reaction
            dfs_traverse(child, next_depth)

    dfs_traverse(route)
    return found_late_amide, findings_json
