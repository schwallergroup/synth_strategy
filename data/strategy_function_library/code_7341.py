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


CROSS_COUPLING_REACTIONS = [
    "Suzuki",
    "Heck",
    "Negishi coupling",
    "Stille reaction",
    "Sonogashira",
    "Buchwald-Hartwig",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the route employs a late-stage cross-coupling reaction in the final three steps.
    This check is performed by testing for a defined set of named reactions, including Suzuki, 
    Heck, Stille, Sonogashira, Buchwald-Hartwig, and others, as specified in the 
    `CROSS_COUPLING_REACTIONS` list.
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

    has_late_cross_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_cross_coupling, findings_json

        if node["type"] == "reaction" and depth <= 2:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for various cross-coupling reactions
                for name in CROSS_COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        has_late_cross_coupling = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break # Found one, no need to check others for this reaction node

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    
    if has_late_cross_coupling:
        # Add the structural constraint if a late-stage cross-coupling was found
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cross_coupling_reaction",
                "position": "last_three_stages"
            }
        })

    return has_late_cross_coupling, findings_json
