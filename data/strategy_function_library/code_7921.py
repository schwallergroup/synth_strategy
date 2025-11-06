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
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (depth <= 1) amide formation by checking if a reaction step matches a predefined list of named amide formation reactions. See the AMIDE_FORMATION_REACTIONS list for all checked reaction types.
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

    late_stage_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_found, findings_json

        if late_stage_amide_found:
            return

        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check if the reaction is a known amide formation type
                for reaction_type in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        late_stage_amide_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        # Add structural constraint if found at late stage
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "amide_formation",
                                "position": "late_stage",
                                "qualifier": "Any reaction from the AMIDE_FORMATION_REACTIONS list at depth <= 1."
                            }
                        })
                        return  # Strategy found, stop processing this branch

        # Traverse children with increased depth
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return late_stage_amide_found, findings_json
