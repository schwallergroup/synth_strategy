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


LATE_STAGE_AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with ammonia to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a late-stage reaction (depth <= 1) is an amide formation by checking against a specific list of named reactions.
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

    late_stage_amide = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_amide, findings_json

        if node["type"] == "reaction" and current_depth <= 1:  # Late stage (depth 0 or 1)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific amide formation reaction types
                is_amide_formation = False
                for rxn in LATE_STAGE_AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_amide_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if is_amide_formation:
                    print(f"Found late-stage amide formation at depth {current_depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    late_stage_amide = True
                    # Add the structural constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "targets": [
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Acyl chloride with ammonia to amide",
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                "Acyl chloride with secondary amine to amide",
                                "Carboxylic acid with primary amine to amide",
                                "Ester with ammonia to amide",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide",
                                "Schotten-Baumann_amide",
                                "Acylation of primary amines",
                                "Acylation of secondary amines",
                                "Carboxylic acid to amide conversion"
                            ],
                            "position": "late_stage",
                            "condition": "depth <= 1"
                        }
                    })

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (which are chemicals)
                dfs_traverse(child, current_depth)
            else:
                # If current node is a chemical, depth increases for children (which are reactions)
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation strategy detected: {late_stage_amide}")
    return late_stage_amide, findings_json
