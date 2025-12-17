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
    This function detects a late-stage SNAr reaction strategy where a nucleophile
    (typically amine) replaces a halogen on a heterocycle.
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

    # Initialize tracking variables
    snar_reaction = False
    snar_depth = float("inf")  # To track how late the SNAr occurs
    max_depth = 0  # To track the total synthesis depth

    def dfs_traverse(node, depth=0):
        nonlocal snar_reaction, snar_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for SNAr reactions using various patterns
                is_snar = False

                # Check for heteroaromatic nucleophilic substitution
                if checker.check_reaction("heteroaromatic_nuc_sub", rsmi):
                    is_snar = True
                    findings_json["atomic_checks"]["named_reactions"].append("heteroaromatic_nuc_sub")
                    print(f"Detected heteroaromatic nucleophilic substitution at depth {depth}")

                # Check for nucleophilic substitution with nitro group activation
                elif checker.check_reaction(
                    "nucl_sub_aromatic_ortho_nitro", rsmi
                ):
                    is_snar = True
                    findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_ortho_nitro")
                    print(f"Detected nitro-activated SNAr at depth {depth}")
                elif checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi):
                    is_snar = True
                    findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_para_nitro")
                    print(f"Detected nitro-activated SNAr at depth {depth}")

                if is_snar:
                    snar_reaction = True
                    snar_depth = min(snar_depth, depth)  # Lower depth means later in synthesis

        # Process children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases for reaction child
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if SNAr occurred in the last third of the synthesis
    is_late_stage = False
    if snar_reaction and max_depth > 0:
        # Lower depth values are later in the synthesis (closer to final product)
        # Consider it late-stage if it's in the first third of the synthesis depth
        is_late_stage = snar_depth <= max_depth / 3
        if is_late_stage:
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "targets": [
                        "heteroaromatic_nuc_sub",
                        "nucl_sub_aromatic_ortho_nitro",
                        "nucl_sub_aromatic_para_nitro"
                    ],
                    "position": "late_stage",
                    "definition": "The reaction must occur in the final third of the synthesis, defined as reaction_depth <= max_synthesis_depth / 3, where depth 0 is the final product."
                }
            })
        print(f"SNAr reaction found at depth {snar_depth} out of maximum depth {max_depth}")
    else:
        print(f"No SNAr reaction found in the synthesis route")

    print(f"Detected late-stage SNAr strategy: {is_late_stage} (depth: {snar_depth}/{max_depth})")
    return is_late_stage, findings_json
