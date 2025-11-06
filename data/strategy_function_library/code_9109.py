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


# Refactoring for Enumeration: Isolate the list of SNAr reaction types.
SNAR_REACTION_TYPES = [
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "heteroaromatic_nuc_sub",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving sequential SNAr reactions on a pyrazole scaffold,
    culminating in a late-stage ether formation. The specific SNAr reaction
    types are defined in the SNAR_REACTION_TYPES list.
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

    # Track key features
    has_pyrazole = False
    snAr_reactions_count = 0
    late_stage_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole, snAr_reactions_count, late_stage_ether_formation, findings_json

        if node["type"] == "mol":
            # Check for pyrazole scaffold
            if node["smiles"]:
                mol_smiles = node["smiles"]
                if checker.check_ring("pyrazole", mol_smiles):
                    has_pyrazole = True
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")

        elif node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for SNAr reactions using reaction checkers
                is_snar = False
                for reaction_type in SNAR_REACTION_TYPES:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_snar = True
                        snAr_reactions_count += 1
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)

                if is_snar:
                    # Check for late-stage ether formation (depth <= 1)
                    # This check is now robust, ensuring the ether is FORMED in this step.
                    if depth <= 1 and checker.check_fg_formation("Ether", rsmi):
                        late_stage_ether_formation = True
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")
                        if "ether_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ether_formation")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases only when moving from chemical to reaction
            next_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if this strategy is present
    strategy_present = has_pyrazole and snAr_reactions_count >= 2 and late_stage_ether_formation

    # Record structural constraints if met
    if has_pyrazole:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "pyrazole"
                ]
            }
        })
    if snAr_reactions_count >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "heteroaromatic_nuc_sub"
                ],
                "operator": ">=",
                "value": 2
            }
        })
    if late_stage_ether_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ether_formation",
                "position": "depth <= 1"
            }
        })

    return strategy_present, findings_json