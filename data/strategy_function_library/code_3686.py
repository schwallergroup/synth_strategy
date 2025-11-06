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


SATURATED_N_HETEROCYCLES_OF_INTEREST = [
    "azetidine",
    "pyrrolidine",
    "piperidine",
    "azepane",
    "morpholine",
    "piperazine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage formation of a specific saturated N-heterocycle via reductive amination. The list of target heterocycles includes azetidine, pyrrolidine, piperidine, azepane, morpholine, and piperazine.
    """
    print("Starting late_stage_reductive_amination function")
    reductive_amination_found = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found, findings_json
        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")

                try:
                    # This single condition correctly checks for the formation of a target
                    # heterocycle via a reductive amination reaction.
                    is_reductive_amination = False
                    if checker.check_reaction("Reductive amination with aldehyde", rsmi):
                        is_reductive_amination = True
                        if "Reductive amination with aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reductive amination with aldehyde")
                    if checker.check_reaction("Reductive amination with ketone", rsmi):
                        is_reductive_amination = True
                        if "Reductive amination with ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Reductive amination with ketone")

                    ring_formed = False
                    rings_found_in_product = []
                    for ring in SATURATED_N_HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, product):
                            ring_formed = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            rings_found_in_product.append(ring)

                    ring_in_reactants = False
                    for r in reactants:
                        for ring in SATURATED_N_HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring, r):
                                ring_in_reactants = True
                                break
                        if ring_in_reactants: break

                    if is_reductive_amination and ring_formed and not ring_in_reactants:
                        print(
                            f"Found late-stage reductive amination with cyclic amine formation at depth {depth}"
                        )
                        reductive_amination_found = True
                        # Add structural constraints
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "Reductive amination with aldehyde",
                                    "Reductive amination with ketone",
                                    "ring_formation"
                                ],
                                "note": "A single reaction must be a ring_formation and one of the two specified reductive amination types."
                            }
                        })
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "ring_formation",
                                "position": "late_stage"
                            }
                        })

                except Exception as e:
                    print(f"Error processing reaction SMILES: {rsmi}")
                    print(f"Error details: {str(e)}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node['type'] != 'reaction' else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    print(f"Reductive amination found: {reductive_amination_found}")

    return reductive_amination_found, findings_json
