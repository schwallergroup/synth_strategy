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


SNAR_REACTIONS_FOR_ETHER_FORMATION = [
    "Nucleophilic aromatic substitution",
    "heteroaromatic_nuc_sub",
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage alkoxy substitution,
    identified by specific named reactions. Late-stage is defined as occurring
    within the final two-thirds of the synthetic sequence. This check is
    triggered by Williamson Ether Synthesis or specific Nucleophilic Aromatic
    Substitution reactions that form an ether, including: Nucleophilic aromatic
    substitution, heteroaromatic_nuc_sub, nucl_sub_aromatic_ortho_nitro, and
    nucl_sub_aromatic_para_nitro.
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

    # Flag to track if we found the alkoxy substitution in the final step
    found_alkoxy_substitution = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    get_max_depth(route)

    # Define what "late-stage" means (e.g., within the last 2/3 of the synthesis depth)
    late_stage_threshold = max_depth * 2 // 3

    def dfs_traverse(node, depth=0):
        nonlocal found_alkoxy_substitution, findings_json

        # Check if this is a late-stage reaction
        is_late_stage = depth <= late_stage_threshold

        if node["type"] == "reaction" and is_late_stage:
            # Get reactants and product
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                parts = rsmi.split(">")
                if len(parts) < 3:
                    return

                reactants_smiles = parts[0].split(".")
                product_smiles = parts[-1]

                # Check if this is a Williamson ether synthesis or similar alkoxy substitution
                if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                    print(f"Found Williamson Ether Synthesis at depth {depth}")
                    found_alkoxy_substitution = True
                    findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Williamson Ether Synthesis",
                            "position": "first_two_thirds_of_retrosynthesis"
                        }
                    })
                    return

                # Check for other types of alkoxy substitution
                ether_formed = checker.check_fg("Ether", product_smiles)
                if ether_formed:
                    if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ether")

                # Check if any reactant has a halide that could be replaced
                has_aromatic_halide = False
                for r in reactants_smiles:
                    if r and checker.check_fg("Aromatic halide", r):
                        has_aromatic_halide = True
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        break

                if has_aromatic_halide and ether_formed:
                    # Check for nucleophilic aromatic substitution reactions
                    for r_name in SNAR_REACTIONS_FOR_ETHER_FORMATION:
                        if checker.check_reaction(r_name, rsmi):
                            print(f"Found nucleophilic aromatic substitution at depth {depth}")
                            found_alkoxy_substitution = True
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)

                            # Add the corresponding structural constraint
                            if r_name == "Nucleophilic aromatic substitution":
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "targets": [
                                            "Nucleophilic aromatic substitution",
                                            "Ether",
                                            "Aromatic halide"
                                        ],
                                        "scope": "reaction_step",
                                        "position": "first_two_thirds_of_retrosynthesis"
                                    }
                                })
                            elif r_name == "heteroaromatic_nuc_sub":
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "targets": [
                                            "heteroaromatic_nuc_sub",
                                            "Ether",
                                            "Aromatic halide"
                                        ],
                                        "scope": "reaction_step",
                                        "position": "first_two_thirds_of_retrosynthesis"
                                    }
                                })
                            elif r_name == "nucl_sub_aromatic_ortho_nitro":
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "targets": [
                                            "nucl_sub_aromatic_ortho_nitro",
                                            "Ether",
                                            "Aromatic halide"
                                        ],
                                        "scope": "reaction_step",
                                        "position": "first_two_thirds_of_retrosynthesis"
                                    }
                                })
                            elif r_name == "nucl_sub_aromatic_para_nitro":
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "targets": [
                                            "nucl_sub_aromatic_para_nitro",
                                            "Ether",
                                            "Aromatic halide"
                                        ],
                                        "scope": "reaction_step",
                                        "position": "first_two_thirds_of_retrosynthesis"
                                    }
                                })
                            return
            except Exception as e:
                print(f"Error in alkoxy substitution detection: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_alkoxy_substitution, findings_json
