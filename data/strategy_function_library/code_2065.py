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


AMIDE_FORMATION_REACTIONS = [
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Ester with ammonia to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final step in the synthesis is an amide bond formation, identified by matching against a predefined list of relevant named reactions and reaction classes.
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

    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation, findings_json

        # Check if this is a reaction node at depth 1 (final synthetic step)
        if node["type"] == "reaction" and depth == 1:
            # Extract reaction SMILES
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check for amide formation using reaction checkers
            for reaction_type in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    found_amide_formation = True
                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    # Add the structural constraint if an amide formation is found at the final stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "targets": AMIDE_FORMATION_REACTIONS,
                            "position": "last_stage"
                        }
                    })
                    return

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # Only increase depth when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            # Optimization: if we already found the feature, no need to traverse further
            if found_amide_formation:
                break
            dfs_traverse(child, next_depth)

    # Check if route is valid
    if not route or "type" not in route:
        return False, findings_json

    # Start traversal
    dfs_traverse(route)

    return found_amide_formation, findings_json
