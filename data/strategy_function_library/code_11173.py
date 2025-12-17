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


HALOGENATION_REACTIONS = [
    "Aromatic fluorination",
    "Aromatic chlorination",
    "Aromatic bromination",
    "Aromatic iodination",
    "Chlorination",
    "Fluorination",
    "Iodination",
    "Bromination",
    "Wohl-Ziegler bromination benzyl primary",
    "Wohl-Ziegler bromination benzyl secondary",
    "Wohl-Ziegler bromination benzyl tertiary",
    "Wohl-Ziegler bromination allyl primary",
    "Wohl-Ziegler bromination allyl secondary",
    "Wohl-Ziegler bromination allyl tertiary",
    "Wohl-Ziegler bromination carbonyl primary",
    "Wohl-Ziegler bromination carbonyl secondary",
    "Wohl-Ziegler bromination carbonyl tertiary",
]

DEHALOGENATION_REACTIONS = ["Aromatic dehalogenation", "Dehalogenation"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a strategy involving halogenation followed by dehalogenation.
    In forward synthesis, halogenation occurs first, then dehalogenation.
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

    halogenation_reactions = []
    dehalogenation_reactions = []

    def dfs_traverse(node, current_depth=0):
        nonlocal halogenation_reactions, dehalogenation_reactions, findings_json
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for halogenation reactions using checker
                for halogenation_type in HALOGENATION_REACTIONS:
                    if checker.check_reaction(halogenation_type, rsmi):
                        halogenation_reactions.append((current_depth, rsmi))
                        if halogenation_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(halogenation_type)
                        break

                # Check for dehalogenation reactions using checker
                for dehalogenation_type in DEHALOGENATION_REACTIONS:
                    if checker.check_reaction(dehalogenation_type, rsmi):
                        dehalogenation_reactions.append((current_depth, rsmi))
                        if dehalogenation_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(dehalogenation_type)
                        break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, current_depth)
            else:  # node['type'] is 'chemical'
                dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)

    has_sequence = False
    if halogenation_reactions and dehalogenation_reactions:
        # Add co-occurrence constraint if both types of reactions are found
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_halogenation_reaction",
                    "any_dehalogenation_reaction"
                ],
                "description": "The route must contain at least one reaction from the halogenation group and at least one from the dehalogenation group."
            }
        })

        halogenation_reactions.sort()
        dehalogenation_reactions.sort()

        for hal_depth, _ in halogenation_reactions:
            for dehal_depth, _ in dehalogenation_reactions:
                if hal_depth > dehal_depth:
                    has_sequence = True
                    # Add sequence constraint if the condition is met
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "before": "any_halogenation_reaction",
                            "after": "any_dehalogenation_reaction",
                            "direction": "forward_synthesis",
                            "description": "In the forward synthesis direction, a halogenation reaction must occur before a dehalogenation reaction."
                        }
                    })
                    break
            if has_sequence:
                break

    return has_sequence, findings_json
