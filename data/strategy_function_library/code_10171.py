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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis includes a temporary epoxide formation and opening sequence.
    This is characterized by formation of an epoxide followed by its opening.
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

    epoxide_formation_depth = None
    epoxide_opening_depth = None
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_formation_depth, epoxide_opening_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for epoxide formation
                reactants_has_epoxide = checker.check_ring("oxirane", reactants_part)
                product_has_epoxide = checker.check_ring("oxirane", product_part)

                # Epoxide formation: no epoxide in reactants, but epoxide in product
                if not reactants_has_epoxide and product_has_epoxide:
                    if epoxide_formation_depth is None or depth < epoxide_formation_depth:
                        epoxide_formation_depth = depth
                        if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Epoxide opening: epoxide in reactants, but no epoxide in product
                if reactants_has_epoxide and not product_has_epoxide:
                    if epoxide_opening_depth is None or depth < epoxide_opening_depth:
                        epoxide_opening_depth = depth
                        if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
            except Exception:
                # Silently ignore errors in malformed reaction data
                pass

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Check if we found both formation and opening, and opening comes after formation in the synthesis
    # (remember lower depth = later in synthesis)
    if epoxide_formation_depth is not None and epoxide_opening_depth is not None:
        if epoxide_opening_depth < epoxide_formation_depth:
            result = True
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "formation_of_oxirane",
                        "opening_of_oxirane"
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "formation_of_oxirane",
                    "after": "opening_of_oxirane"
                }
            })

    return result, findings_json