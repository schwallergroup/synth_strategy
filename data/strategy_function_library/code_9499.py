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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving late-stage amide formation. It specifically checks for the following reaction types: Acylation of Nitrogen Nucleophiles by Carboxylic Acids, and Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N.
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

    found_late_amide_formation = False

    # Calculate depths if not already present
    def calculate_depths(node, current_depth=0):
        node["depth"] = current_depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                calculate_depths(child, current_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                calculate_depths(child, current_depth + 1)

    calculate_depths(route)

    def dfs_traverse(node):
        nonlocal found_late_amide_formation, findings_json

        if node["type"] == "reaction" and node.get("depth", 0) <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if this is an amide formation reaction using a minimal, robust set of checkers
                is_amide_reaction = False
                for reaction_name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amide_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if is_amide_reaction:
                    found_late_amide_formation = True
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if found_late_amide_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N"
                ],
                "position": "depth <= 2"
            }
        })

    return found_late_amide_formation, findings_json
