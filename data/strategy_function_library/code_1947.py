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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acyl chloride with ammonia to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage amide coupling strategy.
    This is confirmed by checking if the final reaction (depth=1) matches one of the predefined named reactions from the `AMIDE_COUPLING_REACTIONS` list.
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

    amide_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_detected, findings_json

        if node["type"] == "reaction":
            # A late-stage reaction is defined as the final step, which is at depth=1
            if depth == 1:
                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    # Check for specific, robustly defined amide coupling reactions
                    for reaction_name in AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            amide_coupling_detected = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            # Add structural constraint if detected at depth 1
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "targets": AMIDE_COUPLING_REACTIONS,
                                    "position": "last_stage"
                                }
                            })
                            break
                except (KeyError, IndexError):
                    # Silently ignore reactions without a valid rsmi
                    pass

        # Traverse children
        for child in node.get("children", []):
            if amide_coupling_detected:
                # Optimization: stop traversing if we've already found the strategy
                return
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return amide_coupling_detected, findings_json
