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


# This list contains general, non-redundant checker names for amide formation.
AMIDE_FORMATION_REACTIONS = [
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves multiple amide bond formations
    at different stages of the synthesis by checking for the following specific reaction types: 'Acylation of primary amines', 'Acylation of secondary amines'.
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

    # Track depths where amide formations occur
    amide_formation_depths = set()
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check for amide formation reactions using a curated list of checkers
            for rxn_name in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"Found amide formation at depth {current_depth}")
                    amide_formation_depths.add(current_depth)
                    if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                    break # Assuming only one match per reaction node is sufficient for recording

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found multiple amide formations at different depths
    if len(amide_formation_depths) >= 2:
        print(f"Detected multiple amide formations at depths: {amide_formation_depths}")
        result = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "amide_formation",
                "operator": ">=",
                "value": 2,
                "count_on": "distinct_stages"
            }
        })

    return result, findings_json
