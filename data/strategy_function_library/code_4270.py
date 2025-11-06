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
    Detects if the synthetic route employs a late-stage reductive amination
    with morpholine as the final step.
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

    # Track if we found the pattern
    found_reductive_amination = False
    reductive_amination_found_flag = False
    morpholine_found_flag = False
    late_stage_reductive_amination_flag = False

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination, reductive_amination_found_flag, morpholine_found_flag, late_stage_reductive_amination_flag, findings_json

        # Check reaction nodes at the late stage (depth <= 1, typically the final step is depth=1)
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if the reaction is a reductive amination using trusted checkers
                is_reductive_amination = False
                ra_reactions = [
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol"
                ]
                for ra_name in ra_reactions:
                    if checker.check_reaction(ra_name, rsmi):
                        is_reductive_amination = True
                        reductive_amination_found_flag = True
                        if ra_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(ra_name)
                        break

                if is_reductive_amination:
                    if depth <= 1:
                        late_stage_reductive_amination_flag = True
                        # Add positional constraint if it's a late-stage reductive amination
                        if {"type": "positional", "details": {"target": "reductive_amination", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "reductive_amination", "position": "late_stage"}})

                    # If it is a reductive amination, check if morpholine is a reactant
                    has_morpholine_reactant = False
                    for reactant in reactants:
                        if checker.check_ring("morpholine", reactant):
                            has_morpholine_reactant = True
                            morpholine_found_flag = True
                            if "morpholine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("morpholine")
                            break
                    
                    if has_morpholine_reactant and depth <= 1:
                        found_reductive_amination = True

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check for co-occurrence constraint after traversal
    if reductive_amination_found_flag and morpholine_found_flag:
        if {"type": "co-occurrence", "details": {"targets": ["reductive_amination", "morpholine"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["reductive_amination", "morpholine"]}})

    return found_reductive_amination, findings_json
