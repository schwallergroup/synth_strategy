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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Ester with secondary amine to amide",
    "Acyl chloride with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Ester with ammonia to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage amide bond formation as the final step by checking if the reaction matches a pre-defined list of named reactions, such as 'Carboxylic acid with primary amine to amide' or 'Schotten-Baumann_amide'.
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

    late_stage_amide_formation = False
    root_is_molecule = route["type"] == "mol"

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation, findings_json

        # Check if this is a reaction at the appropriate depth for late-stage
        is_late_stage_reaction = node["type"] == "reaction" and (
            (root_is_molecule and depth == 1) or (not root_is_molecule and depth == 0)
        )

        if is_late_stage_reaction:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for specific amide formation reaction types
                found_reaction_names = []
                for reaction_name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        found_reaction_names.append(reaction_name)
                
                if found_reaction_names:
                    late_stage_amide_formation = True
                    for r_name in found_reaction_names:
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

            except Exception:
                # Silently ignore errors in reaction analysis
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    if late_stage_amide_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_bond_formation",
                "position": "last_stage"
            }
        })

    return late_stage_amide_formation, findings_json
