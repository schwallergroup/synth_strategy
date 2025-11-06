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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves early formation of an indole core
    followed by subsequent functionalization.
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

    indole_formed = False
    indole_formation_depth = -1.0
    max_depth = 0
    indole_in_starting_material = False

    def dfs_traverse(node, depth=0):
        nonlocal indole_formed, indole_formation_depth, max_depth, indole_in_starting_material, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                product_has_indole = checker.check_ring("indole", product)
                reactants_have_indole = any(checker.check_ring("indole", r) for r in reactants)

                if product_has_indole:
                    if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("indole")

                if product_has_indole and not reactants_have_indole:
                    indole_formed = True
                    # To find the earliest formation, we need the maximum depth
                    indole_formation_depth = max(indole_formation_depth, depth)
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception:
                # Silently ignore errors in reaction analysis
                pass

        elif node["type"] == "mol" and node.get("smiles"):
            if checker.check_ring("indole", node["smiles"]):
                if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("indole")
                if node.get("in_stock", False):
                    indole_in_starting_material = True

        for child in node.get("children", []):
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    early_formation = False
    if indole_formed and not indole_in_starting_material:
        # Early formation means indole is formed in the first third of the synthesis.
        # Larger depth means earlier in the synthesis.
        is_early = indole_formation_depth >= (2 * max_depth / 3)
        early_formation = is_early

    if indole_formed:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ring_formation"
                ],
                "scope": "indole"
            }
        })
    if not indole_in_starting_material:
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "indole",
                "scope": "starting_material"
            }
        })
    if early_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "first_third"
            }
        })

    return early_formation, findings_json
