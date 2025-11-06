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


EARLY_STAGE_POLYCYCLES = [
    "naphthalene",
    "anthracene",
    "phenanthrene",
    "acridine",
    "carbazole",
    "fluorene",
    "indole",
    "quinoline",
    "isoquinoline",
    "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear synthesis that forms a specific polycyclic aromatic system early and maintains it. The strategy is confirmed if the product of the earliest reaction contains a ring from a predefined list (e.g., naphthalene, indole, quinoline) and that same ring system is present in all subsequent products. The synthesis is also checked for linearity (single product per step).
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
    result = False

    # Track reactions and their depths
    reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                reactions.append((node["metadata"]["mapped_reaction_smiles"], depth))
                # Add reaction to atomic checks
                if "reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("reaction")

        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Sort reactions by depth (highest depth = earliest in synthesis)
    reactions.sort(key=lambda x: x[1], reverse=True)

    if len(reactions) < 2:
        print("Not enough reactions in the route (minimum 2 required)")
        return False, findings_json # Return early with findings_json
    else:
        # Structural constraint: at least 2 reactions
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction",
                "operator": ">=",
                "value": 2
            }
        })

    # Check if earliest reaction introduces a polycyclic aromatic system
    earliest_rsmi = reactions[0][0]
    product = earliest_rsmi.split(">")[-1]

    # Check if product contains any polycyclic aromatic system
    has_polycyclic = False
    detected_rings = []

    for ring in EARLY_STAGE_POLYCYCLES:
        if checker.check_ring(ring, product):
            has_polycyclic = True
            detected_rings.append(ring)
            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append(ring)

    if not has_polycyclic:
        print(f"No polycyclic aromatic system found in the earliest reaction product: {product}")
        return False, findings_json # Return early with findings_json
    else:
        # Structural constraint: polycyclic system in first stage product
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": [
                    "naphthalene",
                    "anthracene",
                    "phenanthrene",
                    "acridine",
                    "carbazole",
                    "fluorene",
                    "indole",
                    "quinoline",
                    "isoquinoline",
                    "purine"
                ],
                "entity_type": "ring_system",
                "position": "first_stage",
                "location": "product"
            }
        })

    print(f"Found polycyclic aromatic system(s) in earliest reaction: {', '.join(detected_rings)}")

    # Check if this polycyclic system is maintained in subsequent steps
    polycycle_maintained = True
    for rsmi, depth in reactions[1:]:
        next_product = rsmi.split(">")[-1]

        # Check if any of the detected rings are still present
        rings_present = False
        for ring in detected_rings:
            if checker.check_ring(ring, next_product):
                rings_present = True
                break

        if not rings_present:
            print(f"Polycyclic aromatic system lost at depth {depth}, product: {next_product}")
            polycycle_maintained = False
            # Add loss of early stage polycycle to atomic checks
            if "loss_of_early_stage_polycycle" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("loss_of_early_stage_polycycle")
            break

    if not polycycle_maintained:
        return False, findings_json # Return early with findings_json
    else:
        # Structural constraint: negation of loss of early stage polycycle
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "loss_of_early_stage_polycycle",
                "scope": "any_subsequent_step"
            }
        })

    # Check if synthesis is linear (each reaction has only one product)
    is_linear = True
    for rsmi, depth in reactions:
        products = rsmi.split(">")[-1].split(".")
        if len(products) > 1:
            print(f"Non-linear synthesis detected at depth {depth}: multiple products found")
            is_linear = False
            # Add multi-product reaction to atomic checks
            if "multi-product_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                findings_json["atomic_checks"]["named_reactions"].append("multi-product_reaction")
            break

    if not is_linear:
        return False, findings_json # Return early with findings_json
    else:
        # Structural constraint: negation of multi-product reaction
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "multi-product_reaction",
                "scope": "any_step"
            }
        })

    print("Found linear synthesis maintaining polycyclic aromatic system")
    result = True
    return result, findings_json
