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

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if the route contains late-stage epoxide opening.
    Late stage means it happens at a low depth in the retrosynthetic tree.
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

    epoxide_opening_reactions = []
    result = False

    def dfs_traverse(node, depth=0, path=[]):
        nonlocal epoxide_opening_reactions, findings_json
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rxn_smiles = node["metadata"]["mapped_reaction_smiles"]

                # Check for epoxide opening by functional group changes
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check if any reactant contains an epoxide
                epoxide_in_reactants = False
                for reactant in reactants:
                    if checker.check_ring("oxirane", reactant):
                        epoxide_in_reactants = True
                        if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                        break

                # Check if the product doesn't contain an epoxide (it was opened)
                epoxide_in_product = checker.check_ring("oxirane", product)

                if epoxide_in_reactants and not epoxide_in_product:
                    epoxide_opening_reactions.append((depth, rxn_smiles))
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

        # Continue traversal
        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth, path + [node])

    dfs_traverse(route)

    # Check if we have epoxide opening reactions
    if not epoxide_opening_reactions:
        return False, findings_json

    # Check if any epoxide opening happens at late stage (low depth)
    # Let's define late stage as depth <= 3 to be more flexible
    late_stage_openings = [rxn for depth, rxn in epoxide_opening_reactions if depth <= 3]

    if late_stage_openings:
        result = True
        # Add the structural constraint if late-stage opening is found
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_destruction",
                "context": "oxirane",
                "position": "late_stage",
                "condition": "depth <= 3"
            }
        })

    return result, findings_json